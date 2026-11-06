!===============================================================================
! CUBE_FoF: Friends-of-Friends (FOF) halo finder with Adaptive Bound Test
!
! 功能说明:
!   对给定的粒子相空间与质量执行两遍空间 CRS + 13 半邻域 + 并查集 (DSU) 找链，
!   并按需调用多分辨率自适应 PM + C2 短程自洽去束缚批处理求解器。
!   支持单/多区域独立找链、变质量权重、集中式并行核心调度与无进程崩溃异常捕获。
!===============================================================================

subroutine fof_fortran(xp_halo, vp_halo, mp_halo, np_halo, nhalo, L_Box, &
   linking_length, np_halo_min, bound_check, verbose, softening, &
   a_val, b_val, s_15_val, n_threads, total_np, &
   in_halo, halo_id, ang_mom, ierr, ng_val, ng_global_val)
   use bound_test, only: bound_set_context, bound_check_batch, bound_release, &
      xv, mp, pid, bound_member_table, &
      bound_member_offset, bound_member_count, fofcore
   use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads, &
      omp_set_max_active_levels, omp_in_parallel, omp_set_dynamic
   use iso_c_binding, only: c_int32_t
   implicit none

   !---------------------------------------------------------------------------
   ! 输入参数声明
   !---------------------------------------------------------------------------
   integer, intent(in) :: nhalo                          ! halo/区域数量
   integer, intent(in) :: total_np                       ! 总粒子数
   integer, intent(in) :: np_halo(nhalo)                 ! 每个区域的粒子数
   real(8), intent(in) :: xp_halo(3, total_np)           ! 粒子坐标 (3, total_np)
   real(8), intent(in) :: vp_halo(3, total_np)           ! 粒子速度 (3, total_np)
   real(8), intent(in) :: mp_halo(total_np)              ! 粒子质量权重 (total_np)
   real(8), intent(in) :: L_Box                          ! 模拟盒边长
   real(8), intent(in) :: linking_length                 ! FOF链接长度
   integer, intent(in) :: np_halo_min                    ! 最小粒子数阈值
   integer, intent(in) :: bound_check                    ! 是否束缚检查 (0=False, 1=True)
   integer, intent(in) :: verbose                        ! 详细输出级别
   real(8), intent(in) :: softening                      ! 软化长度
   real(8), intent(in) :: a_val                          ! 尺度因子 a
   real(8), intent(in) :: b_val                          ! 哈勃相对流系数 B
   real(8), intent(in) :: s_15_val                       ! 15 Mpc 对应的格点跨度
   integer, intent(in) :: n_threads                      ! 允许的最大线程数
   integer(8), optional, intent(in) :: ng_val            ! 单节点网格维度 ng (默认 256)
   real(8), optional, intent(in) :: ng_global_val        ! 全局网格维度 ng_global (默认 512.0)

   !---------------------------------------------------------------------------
   ! 输出参数声明
   !---------------------------------------------------------------------------
   integer, intent(out) :: in_halo(total_np)             ! 输出：粒子是否在有效halo中 (0/1)
   integer, intent(out) :: halo_id(total_np)             ! 输出：粒子所属有效halo编号 (0表示不属于)
   real(8), intent(out) :: ang_mom(3, nhalo)             ! 输出：每个区域累加的角动量
   integer, intent(out) :: ierr                          ! 输出：错误码 (0=正常)

   !f2py intent(in) :: nhalo, total_np, np_halo, xp_halo, vp_halo, mp_halo
   !f2py intent(in) :: L_Box, linking_length, np_halo_min, bound_check, verbose
   !f2py intent(in) :: softening, a_val, b_val, s_15_val, n_threads
   !f2py intent(in), optional :: ng_val, ng_global_val
   !f2py intent(out) :: in_halo, halo_id, ang_mom, ierr
   !f2py depend(nhalo) :: np_halo, ang_mom
   !f2py depend(total_np) :: xp_halo, vp_halo, mp_halo, in_halo, halo_id

   !---------------------------------------------------------------------------
   ! C 无锁原子并查集接口
   !---------------------------------------------------------------------------
   interface
      subroutine fof_atomic_union(parent, i, j) bind(c, name="fof_atomic_union")
         import :: c_int32_t
         integer(c_int32_t), intent(inout) :: parent(*)
         integer(c_int32_t), value :: i, j
      end subroutine fof_atomic_union
   end interface

   !---------------------------------------------------------------------------
   ! 13 半邻域模版
   !---------------------------------------------------------------------------
   integer, parameter :: ijk_stencil(3, 13) = reshape([ &
      -1,-1,-1,  0,-1,-1,  1,-1,-1, &
      -1, 0,-1,  0, 0,-1,  1, 0,-1, &
      -1, 1,-1,  0, 1,-1,  1, 1,-1, &
      -1,-1, 0,  0,-1, 0,  1,-1, 0, &
      -1, 0, 0 ], [3, 13])

   !---------------------------------------------------------------------------
   ! 区域候选缓冲派生类型
   !---------------------------------------------------------------------------
   type region_cand_type
      integer :: ncand = 0
      integer, allocatable :: cand_np(:)
      integer, allocatable :: cand_head(:)
      integer, allocatable :: cand_members(:)
      integer, allocatable :: member_ptr(:)
      real(4), allocatable :: xm(:,:)
      real(4), allocatable :: vm(:,:)
   end type region_cand_type

   type(region_cand_type), allocatable :: reg_cands(:)
   integer, allocatable :: reg_start(:), order_reg(:)
   integer :: T, fof_team, fof_nest, N_max, ihalo, reg_k, istart, np, valid_h
   integer :: nh_cand, total_cand_np, h, p, orig_idx, k, res_offset, g, i
   integer :: anchor_idx, ios, alloc_err_reg
   real(8) :: L_Box2, dx_rel(3), dv_rel(3), wi, w_sum, qc(3), vc(3), xp_anchor_cand(3)

   ! 束缚求解批数组
   integer, allocatable :: heads(:), np_cand_list(:), nb(:), status(:), final_members(:)
   integer, allocatable :: cand_region(:)
   real(4), allocatable :: xm_batch(:,:), vm_batch(:,:), xb_batch(:,:), vb_batch(:,:)
   integer(8) :: ng_eff
   real(8) :: ng_global_eff

   !---------------------------------------------------------------------------
   ! 1. 初始化输出与基础参数校验
   !---------------------------------------------------------------------------
   in_halo = 0
   halo_id = 0
   ang_mom = 0.0d0
   ierr = 0

   if (n_threads > 256) then
      ierr = 10
      return
   end if

   if (nhalo <= 0 .or. total_np <= 0) then
      return
   end if

   if (n_threads <= 0) then
      T = min(omp_get_max_threads(), 256)
   else
      T = min(n_threads, 256)
   end if

   N_max = maxval(np_halo)
   call compute_fof_parallel_config(T, nhalo, N_max, fof_team, fof_nest)

   L_Box2 = 0.5d0 * L_Box

   allocate(reg_start(nhalo), order_reg(nhalo), reg_cands(nhalo), stat=ios)
   if (ios /= 0) then
      ierr = 3
      return
   end if
   reg_start(1) = 1
   do i = 1, nhalo - 1
      reg_start(i + 1) = reg_start(i) + np_halo(i)
   end do

   ! 按区域粒子数从大到小生成排序索引，利于外层动态负载均衡
   do i = 1, nhalo
      order_reg(i) = i
   end do
   call sort_regions_descending(nhalo, np_halo, order_reg)

   !---------------------------------------------------------------------------
   ! 2. 跨区域并行找链（两遍空间 CRS + 13 半邻域 + DSU）
   !---------------------------------------------------------------------------
   call omp_set_max_active_levels(2)
   call omp_set_dynamic(.false.)

   !$omp parallel do num_threads(fof_team) schedule(dynamic, 1) if (fof_team > 1) &
   !$omp private(reg_k, ihalo, istart, np, alloc_err_reg)
   do reg_k = 1, nhalo
      ihalo = order_reg(reg_k)
      istart = reg_start(ihalo)
      np = np_halo(ihalo)
      if (np >= np_halo_min .and. np > 0) then
         call process_single_region(istart, np, L_Box, L_Box2, linking_length, &
            np_halo_min, fof_nest, reg_cands(ihalo), alloc_err_reg)
         if (alloc_err_reg /= 0) then
            !$omp atomic write
            ierr = alloc_err_reg
         end if
      end if
   end do
   !$omp end parallel do

   if (ierr /= 0) then
      deallocate(reg_start, order_reg, reg_cands)
      return
   end if

   !---------------------------------------------------------------------------
   ! 3. 串行汇总所有区域的候选团块
   !---------------------------------------------------------------------------
   nh_cand = 0
   total_cand_np = 0
   do ihalo = 1, nhalo
      if (reg_cands(ihalo)%ncand > 0) then
         nh_cand = nh_cand + reg_cands(ihalo)%ncand
         total_cand_np = total_cand_np + size(reg_cands(ihalo)%cand_members)
      end if
   end do

   if (nh_cand == 0) then
      deallocate(reg_start, order_reg, reg_cands)
      return
   end if

   allocate(cand_region(nh_cand), stat=ios)
   if (ios /= 0) then
      ierr = 3
      deallocate(reg_start, order_reg, reg_cands)
      return
   end if

   !---------------------------------------------------------------------------
   ! 4. 束缚检查与写回
   !---------------------------------------------------------------------------
   if (bound_check == 0) then
      ! 无需束缚检查：直接输出全部候选团块
      h = 0
      do ihalo = 1, nhalo
         if (reg_cands(ihalo)%ncand <= 0) cycle
         do g = 1, reg_cands(ihalo)%ncand
            h = h + 1
            cand_region(h) = ihalo
            anchor_idx = reg_cands(ihalo)%cand_members(reg_cands(ihalo)%member_ptr(g))
            xp_anchor_cand = xp_halo(:, anchor_idx)
            ! 标记成员
            do p = reg_cands(ihalo)%member_ptr(g), reg_cands(ihalo)%member_ptr(g + 1) - 1
               orig_idx = reg_cands(ihalo)%cand_members(p)
               in_halo(orig_idx) = 1
               halo_id(orig_idx) = h
            end do

            ! 质量加权质心与均速，累加角动量
            w_sum = 0.0d0
            qc = 0.0d0
            vc = 0.0d0
            do p = reg_cands(ihalo)%member_ptr(g), reg_cands(ihalo)%member_ptr(g + 1) - 1
               orig_idx = reg_cands(ihalo)%cand_members(p)
               wi = mp_halo(orig_idx)
               w_sum = w_sum + wi
               dx_rel = modulo(xp_halo(:, orig_idx) - xp_anchor_cand + L_Box2, L_Box) - L_Box2
               qc = qc + wi * dx_rel
               vc = vc + wi * vp_halo(:, orig_idx)
            end do
            if (w_sum > 0.0d0) then
               qc = qc / w_sum
               vc = vc / w_sum
               do p = reg_cands(ihalo)%member_ptr(g), reg_cands(ihalo)%member_ptr(g + 1) - 1
                  orig_idx = reg_cands(ihalo)%cand_members(p)
                  wi = mp_halo(orig_idx)
                  dx_rel = (modulo(xp_halo(:, orig_idx) - xp_anchor_cand + L_Box2, L_Box) - L_Box2) - qc
                  dv_rel = vp_halo(:, orig_idx) - vc
                  ang_mom(1, ihalo) = ang_mom(1, ihalo) + wi * (dx_rel(2) * dv_rel(3) - dx_rel(3) * dv_rel(2))
                  ang_mom(2, ihalo) = ang_mom(2, ihalo) + wi * (dx_rel(3) * dv_rel(1) - dx_rel(1) * dv_rel(3))
                  ang_mom(3, ihalo) = ang_mom(3, ihalo) + wi * (dx_rel(1) * dv_rel(2) - dx_rel(2) * dv_rel(1))
               end do
            end if
         end do
      end do

   else
      ! bound_check == 1: 自适应自洽束缚检查求解
      ng_eff = 256_8
      if (present(ng_val)) then
         if (ng_val > 0_8) ng_eff = ng_val
      end if
      ng_global_eff = 512.0d0
      if (present(ng_global_val)) then
         if (ng_global_val > 0.0d0) ng_global_eff = ng_global_val
      end if

      call bound_set_context(softening, a_val, b_val, s_15_val, np_halo_min, verbose, T, ng_eff, ng_global_eff)

      allocate(xv(6, total_cand_np), mp(total_cand_np), pid(total_cand_np), stat=ios)
      if (ios == 0) allocate(bound_member_table(total_cand_np), bound_member_offset(nh_cand), bound_member_count(nh_cand), stat=ios)
      if (ios == 0) allocate(heads(nh_cand), np_cand_list(nh_cand), xm_batch(3, nh_cand), vm_batch(3, nh_cand), stat=ios)
      if (ios == 0) allocate(nb(nh_cand), xb_batch(3, nh_cand), vb_batch(3, nh_cand), status(nh_cand), stat=ios)
      if (ios == 0) allocate(final_members(total_cand_np), stat=ios)
      if (ios /= 0) then
         ierr = 3
         if (allocated(xv)) deallocate(xv, mp, pid)
         if (allocated(bound_member_table)) deallocate(bound_member_table, bound_member_offset, bound_member_count)
         if (allocated(heads)) deallocate(heads, np_cand_list, xm_batch, vm_batch)
         if (allocated(nb)) deallocate(nb, xb_batch, vb_batch, status)
         if (allocated(final_members)) deallocate(final_members)
         deallocate(cand_region, reg_start, order_reg, reg_cands)
         return
      end if

      h = 0
      res_offset = 1
      do ihalo = 1, nhalo
         if (reg_cands(ihalo)%ncand <= 0) cycle
         do g = 1, reg_cands(ihalo)%ncand
            h = h + 1
            cand_region(h) = ihalo
            np_cand_list(h) = reg_cands(ihalo)%cand_np(g)
            bound_member_offset(h) = res_offset
            bound_member_count(h) = np_cand_list(h)
            heads(h) = res_offset
            xm_batch(:, h) = reg_cands(ihalo)%xm(:, g)
            vm_batch(:, h) = reg_cands(ihalo)%vm(:, g)
            anchor_idx = reg_cands(ihalo)%cand_members(reg_cands(ihalo)%member_ptr(g))
            xp_anchor_cand = xp_halo(:, anchor_idx)

            do p = 1, np_cand_list(h)
               orig_idx = reg_cands(ihalo)%cand_members(reg_cands(ihalo)%member_ptr(g) + p - 1)
               dx_rel = modulo(xp_halo(:, orig_idx) - xp_anchor_cand + L_Box2, L_Box) - L_Box2
               xv(1:3, res_offset + p - 1) = real(dx_rel, 4)
               xv(4:6, res_offset + p - 1) = real(vp_halo(:, orig_idx), 4)
               mp(res_offset + p - 1) = real(mp_halo(orig_idx), 4)
               pid(res_offset + p - 1) = orig_idx
               bound_member_table(res_offset + p - 1) = res_offset + p - 1
            end do
            res_offset = res_offset + np_cand_list(h)
         end do
      end do

      call bound_check_batch(heads, np_cand_list, xm_batch, vm_batch, nb, xb_batch, vb_batch, status, &
         max_threads=T, final_members=final_members, ierr=ierr)

      if (ierr == 0) then
         valid_h = 0
         res_offset = 1
         do h = 1, nh_cand
            ihalo = cand_region(h)
            ! 仅 status == 0 且留存数 >= np_halo_min 视为有效束缚 Halo
            if (status(h) == 0 .and. nb(h) >= np_halo_min) then
               valid_h = valid_h + 1
               anchor_idx = int(pid(heads(h)))
               xp_anchor_cand = xp_halo(:, anchor_idx)
               ! 计算保留成员的加权质心与均速
               w_sum = 0.0d0
               qc = 0.0d0
               vc = 0.0d0
               do k = 1, nb(h)
                  p = final_members(res_offset + k - 1)
                  orig_idx = int(pid(p))
                  wi = mp_halo(orig_idx)
                  w_sum = w_sum + wi
                  dx_rel = modulo(xp_halo(:, orig_idx) - xp_anchor_cand + L_Box2, L_Box) - L_Box2
                  qc = qc + wi * dx_rel
                  vc = vc + wi * vp_halo(:, orig_idx)
               end do
               if (w_sum > 0.0d0) then
                  qc = qc / w_sum
                  vc = vc / w_sum
                  do k = 1, nb(h)
                     p = final_members(res_offset + k - 1)
                     orig_idx = int(pid(p))
                     wi = mp_halo(orig_idx)
                     in_halo(orig_idx) = 1
                     halo_id(orig_idx) = valid_h
                     dx_rel = (modulo(xp_halo(:, orig_idx) - xp_anchor_cand + L_Box2, L_Box) - L_Box2) - qc
                     dv_rel = vp_halo(:, orig_idx) - vc
                     ang_mom(1, ihalo) = ang_mom(1, ihalo) + wi * (dx_rel(2) * dv_rel(3) - dx_rel(3) * dv_rel(2))
                     ang_mom(2, ihalo) = ang_mom(2, ihalo) + wi * (dx_rel(3) * dv_rel(1) - dx_rel(1) * dv_rel(3))
                     ang_mom(3, ihalo) = ang_mom(3, ihalo) + wi * (dx_rel(1) * dv_rel(2) - dx_rel(2) * dv_rel(1))
                  end do
               end if
            end if
            res_offset = res_offset + np_cand_list(h)
         end do
      end if

      call bound_release
      deallocate(xv, mp, pid, bound_member_table, bound_member_offset, bound_member_count)
      deallocate(heads, np_cand_list, xm_batch, vm_batch, nb, xb_batch, vb_batch, status, final_members)
   end if

   deallocate(cand_region, reg_start, order_reg, reg_cands)

contains

   !---------------------------------------------------------------------------
   ! 集中式 FOF 并行参数计算
   !---------------------------------------------------------------------------
   subroutine compute_fof_parallel_config(T_in, H_in, N_max_in, fof_team_out, fof_nest_out)
      integer, intent(in) :: T_in, H_in, N_max_in
      integer, intent(out) :: fof_team_out, fof_nest_out

      if (T_in >= 2 .and. H_in <= 10 * T_in .and. N_max_in > 10000) then
         fof_nest_out = 2
      else
         fof_nest_out = 1
      end if
      fof_team_out = min(H_in, max(1, T_in / fof_nest_out))
   end subroutine compute_fof_parallel_config

   !---------------------------------------------------------------------------
   ! 区域排序子程序 (按粒子数降序)
   !---------------------------------------------------------------------------
   subroutine sort_regions_descending(n, counts, order)
      integer, intent(in) :: n, counts(n)
      integer, intent(inout) :: order(n)
      integer :: l, r, i, j, key_idx, key_val

      if (n <= 1) return
      l = n / 2 + 1
      r = n
      do
         if (l > 1) then
            l = l - 1
            key_idx = order(l)
         else
            key_idx = order(r)
            order(r) = order(1)
            r = r - 1
            if (r == 1) then
               order(1) = key_idx
               exit
            end if
         end if
         key_val = counts(key_idx)
         i = l
         j = 2 * l
         do while (j <= r)
            if (j < r) then
               if (counts(order(j + 1)) < counts(order(j))) j = j + 1
            end if
            if (key_val <= counts(order(j))) exit
            order(i) = order(j)
            i = j
            j = 2 * i
         end do
         order(i) = key_idx
      end do
   end subroutine sort_regions_descending

   !---------------------------------------------------------------------------
   ! 串行快速并查集 union
   !---------------------------------------------------------------------------
   subroutine dsu_union_serial(parent, i, j)
      integer, intent(inout) :: parent(*)
      integer, intent(in) :: i, j
      integer :: ra, rb, next_node, curr

      ra = i
      do while (parent(ra) > 0)
         ra = parent(ra)
      end do
      curr = i
      do while (curr /= ra)
         next_node = parent(curr)
         parent(curr) = ra
         curr = next_node
      end do

      rb = j
      do while (parent(rb) > 0)
         rb = parent(rb)
      end do
      curr = j
      do while (curr /= rb)
         next_node = parent(curr)
         parent(curr) = rb
         curr = next_node
      end do

      if (ra == rb) return
      ! 大编号挂在小编号下
      if (ra > rb) then
         parent(ra) = rb
      else
         parent(rb) = ra
      end if
   end subroutine dsu_union_serial

   !---------------------------------------------------------------------------
   ! 查找根节点带路径压缩
   !---------------------------------------------------------------------------
   integer function dsu_find_root(parent, i) result(iroot)
      integer, intent(inout) :: parent(*)
      integer, intent(in) :: i
      integer :: inext, ipart

      iroot = i
      do while (parent(iroot) > 0)
         iroot = parent(iroot)
      end do
      ipart = i
      do while (ipart /= iroot)
         inext = parent(ipart)
         parent(ipart) = iroot
         ipart = inext
      end do
   end function dsu_find_root

   !---------------------------------------------------------------------------
   ! 周期最小镜像距离平方 (两点)
   !---------------------------------------------------------------------------
   pure real(8) function dist2_pbc(x1, x2, box_len, box_half)
      real(4), intent(in) :: x1(3), x2(3)
      real(8), intent(in) :: box_len, box_half
      real(8) :: dx, dy, dz
      dx = abs(dble(x1(1) - x2(1))); if (dx > box_half) dx = box_len - dx
      dy = abs(dble(x1(2) - x2(2))); if (dy > box_half) dy = box_len - dy
      dz = abs(dble(x1(3) - x2(3))); if (dz > box_half) dz = box_len - dz
      dist2_pbc = dx*dx + dy*dy + dz*dz
   end function dist2_pbc

   !---------------------------------------------------------------------------
   ! 处理单个输入区域：两遍 CRS + 13 半邻域/周期CRS + DSU 找链
   !---------------------------------------------------------------------------
   subroutine process_single_region(istart_in, np_in, box_len, box_half, &
      r_link, np_min, nest_threads, cand_out, alloc_err)
      integer, intent(in) :: istart_in, np_in, np_min, nest_threads
      real(8), intent(in) :: box_len, box_half, r_link
      type(region_cand_type), intent(out) :: cand_out
      integer, intent(out) :: alloc_err

      real(8) :: x_rel(3, np_in), xp_anchor(3), xp_min(3), xp_max(3), xp_range(3)
      real(8) :: vol, L_mesh, rp2, dx(3), r2, w_sum, wi, qc_d(3), vc_d(3), xp_anchor_cand(3)
      real(8) :: cell_origin(3), cell_size(3)
      integer :: Nx, Ny, Nz, ncell, icell, i, ix, iy, iz, p, p1, p2, d
      integer :: jx, jy, jz, jcell, i_n, root_p, g, cand_cnt, mem_cnt
      integer :: N_dim(3), dz, dy, dx_step, num_nbr, nbr_cells(27), k, ios, anchor_idx
      logical :: periodic_dim(3), is_periodic, duplicate
      integer, allocatable :: cell_ptr(:), cell_cursor(:), cell_cursor_save(:), cand_cursor(:)
      integer, allocatable :: parent_size(:), pid_reg(:)
      real(4), allocatable :: x_crs(:,:), v_crs(:,:), mp_crs(:)
      integer, allocatable :: root_size(:), root_map(:)
      logical, allocatable :: layer_done(:)
      integer :: my_start, step, ithread
      logical :: claimed

      alloc_err = 0
      cand_out%ncand = 0
      rp2 = r_link * r_link
      xp_anchor = xp_halo(:, istart_in)

      ! 周期展开至以首粒子为中心的相对坐标
      do i = 1, np_in
         x_rel(:, i) = modulo(xp_halo(:, istart_in + i - 1) - xp_anchor + box_half, box_len) - box_half
      end do

      xp_min(1) = minval(x_rel(1, :)); xp_max(1) = maxval(x_rel(1, :))
      xp_min(2) = minval(x_rel(2, :)); xp_max(2) = maxval(x_rel(2, :))
      xp_min(3) = minval(x_rel(3, :)); xp_max(3) = maxval(x_rel(3, :))
      xp_range = xp_max - xp_min

      do d = 1, 3
         periodic_dim(d) = (xp_range(d) + 2.0d0 * r_link >= box_len)
      end do
      is_periodic = any(periodic_dim)

      vol = max(xp_range(1), r_link) * max(xp_range(2), r_link) * max(xp_range(3), r_link)
      L_mesh = max(r_link, (vol / dble(max(1, np_in)))**(1.0d0 / 3.0d0))

      do d = 1, 3
         if (periodic_dim(d)) then
            cell_origin(d) = -box_half
            N_dim(d) = max(1, min(int(box_len / r_link), int(box_len / L_mesh)))
            cell_size(d) = box_len / dble(N_dim(d))
         else
            cell_origin(d) = xp_min(d)
            N_dim(d) = max(1, ceiling(xp_range(d) / L_mesh))
            cell_size(d) = L_mesh
         end if
      end do
      Nx = N_dim(1); Ny = N_dim(2); Nz = N_dim(3)

      if (dble(Nx) * dble(Ny) * dble(Nz) > dble(max(262144, 4 * np_in))) then
         L_mesh = (vol / dble(max(1, np_in) * 2))**(1.0d0 / 3.0d0)
         L_mesh = max(r_link, L_mesh)
         do d = 1, 3
            if (periodic_dim(d)) then
               cell_origin(d) = -box_half
               N_dim(d) = max(1, min(int(box_len / r_link), int(box_len / L_mesh)))
               cell_size(d) = box_len / dble(N_dim(d))
            else
               cell_origin(d) = xp_min(d)
               N_dim(d) = max(1, ceiling(xp_range(d) / L_mesh))
               cell_size(d) = L_mesh
            end if
         end do
         Nx = N_dim(1); Ny = N_dim(2); Nz = N_dim(3)
      end if
      ncell = Nx * Ny * Nz

      allocate(cell_ptr(ncell + 1), cell_cursor(ncell), cell_cursor_save(ncell), stat=ios)
      if (ios /= 0) then
         alloc_err = 3
         return
      end if
      cell_cursor = 0

      ! Pass 1: 统计各 cell 粒子数
      do i = 1, np_in
         ix = min(Nx, max(1, floor((x_rel(1, i) - cell_origin(1)) / cell_size(1)) + 1))
         iy = min(Ny, max(1, floor((x_rel(2, i) - cell_origin(2)) / cell_size(2)) + 1))
         iz = min(Nz, max(1, floor((x_rel(3, i) - cell_origin(3)) / cell_size(3)) + 1))
         icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
         cell_cursor(icell) = cell_cursor(icell) + 1
      end do

      cell_ptr(1) = 1
      do icell = 1, ncell
         cell_ptr(icell + 1) = cell_ptr(icell) + cell_cursor(icell)
         cell_cursor(icell) = cell_ptr(icell)
      end do
      cell_cursor_save = cell_cursor

      ! Pass 2: 装配空间 CRS 连续数组
      allocate(x_crs(3, np_in), v_crs(3, np_in), mp_crs(np_in), pid_reg(np_in), stat=ios)
      if (ios /= 0) then
         alloc_err = 3
         deallocate(cell_ptr, cell_cursor, cell_cursor_save)
         return
      end if
      do i = 1, np_in
         ix = min(Nx, max(1, floor((x_rel(1, i) - cell_origin(1)) / cell_size(1)) + 1))
         iy = min(Ny, max(1, floor((x_rel(2, i) - cell_origin(2)) / cell_size(2)) + 1))
         iz = min(Nz, max(1, floor((x_rel(3, i) - cell_origin(3)) / cell_size(3)) + 1))
         icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
         p = cell_cursor(icell)
         x_crs(:, p) = real(x_rel(:, i), 4)
         v_crs(:, p) = real(vp_halo(:, istart_in + i - 1), 4)
         mp_crs(p) = real(mp_halo(istart_in + i - 1), 4)
         pid_reg(p) = istart_in + i - 1
         cell_cursor(icell) = p + 1
      end do

      ! 初始化并查集
      allocate(parent_size(np_in), stat=ios)
      if (ios /= 0) then
         alloc_err = 3
         deallocate(cell_ptr, cell_cursor, cell_cursor_save)
         deallocate(x_crs, v_crs, mp_crs, pid_reg)
         return
      end if
      parent_size = -1

      if (.not. is_periodic) then
         ! 非周期：快速 13 半邻域
         if (nest_threads <= 1) then
            do iz = 1, Nz
               do iy = 1, Ny
                  do ix = 1, Nx
                     icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
                     if (cell_ptr(icell) >= cell_ptr(icell + 1)) cycle
                     ! 本 cell 内上三角
                     do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                        do p2 = p1 + 1, cell_ptr(icell + 1) - 1
                           r2 = sum((x_crs(:, p1) - x_crs(:, p2))**2)
                           if (r2 <= rp2) call dsu_union_serial(parent_size, p1, p2)
                        end do
                     end do
                     ! 13 半邻域
                     do i_n = 1, 13
                        jx = ix + ijk_stencil(1, i_n)
                        jy = iy + ijk_stencil(2, i_n)
                        jz = iz + ijk_stencil(3, i_n)
                        if (jx < 1 .or. jx > Nx .or. jy < 1 .or. jy > Ny .or. jz < 1 .or. jz > Nz) cycle
                        jcell = jx + Nx * (jy - 1) + Nx * Ny * (jz - 1)
                        if (cell_ptr(jcell) >= cell_ptr(jcell + 1)) cycle
                        do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                           do p2 = cell_ptr(jcell), cell_ptr(jcell + 1) - 1
                              r2 = sum((x_crs(:, p1) - x_crs(:, p2))**2)
                              if (r2 <= rp2) call dsu_union_serial(parent_size, p1, p2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do

         else
            ! 并行 z 层环扫抢占 + 原子并查集
            allocate(layer_done(Nz), stat=ios)
            if (ios /= 0) then
               alloc_err = 3
               deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
               deallocate(x_crs, v_crs, mp_crs)
               return
            end if
            layer_done = .false.
            !$omp parallel num_threads(nest_threads) default(shared) &
            !$omp private(ithread, my_start, step, iz, claimed, iy, ix, icell, p1, p2, r2, i_n, jx, jy, jz, jcell)
            ithread = omp_get_thread_num() + 1
            my_start = int(int(Nz, 8) * (ithread - 1) / omp_get_num_threads()) + 1
            do step = 0, Nz - 1
               iz = mod(my_start - 1 + step, Nz) + 1
               !$omp atomic capture
               claimed = layer_done(iz)
               layer_done(iz) = .true.
               !$omp end atomic
               if (claimed) cycle

               do iy = 1, Ny
                  do ix = 1, Nx
                     icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
                     if (cell_ptr(icell) >= cell_ptr(icell + 1)) cycle
                     ! 本 cell 内上三角
                     do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                        do p2 = p1 + 1, cell_ptr(icell + 1) - 1
                           r2 = sum((x_crs(:, p1) - x_crs(:, p2))**2)
                           if (r2 <= rp2) call fof_atomic_union(parent_size, p1, p2)
                        end do
                     end do
                     ! 13 半邻域
                     do i_n = 1, 13
                        jx = ix + ijk_stencil(1, i_n)
                        jy = iy + ijk_stencil(2, i_n)
                        jz = iz + ijk_stencil(3, i_n)
                        if (jx < 1 .or. jx > Nx .or. jy < 1 .or. jy > Ny .or. jz < 1 .or. jz > Nz) cycle
                        jcell = jx + Nx * (jy - 1) + Nx * Ny * (jz - 1)
                        if (cell_ptr(jcell) >= cell_ptr(jcell + 1)) cycle
                        do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                           do p2 = cell_ptr(jcell), cell_ptr(jcell + 1) - 1
                              r2 = sum((x_crs(:, p1) - x_crs(:, p2))**2)
                              if (r2 <= rp2) call fof_atomic_union(parent_size, p1, p2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
            !$omp end parallel
            deallocate(layer_done)
         end if

      else
         ! 周期：支持跨周期卷绕与最小镜像距离 dist2_pbc
         if (nest_threads <= 1) then
            do iz = 1, Nz
               do iy = 1, Ny
                  do ix = 1, Nx
                     icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
                     if (cell_ptr(icell) >= cell_ptr(icell + 1)) cycle
                     ! 本 cell 内上三角
                     do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                        do p2 = p1 + 1, cell_ptr(icell + 1) - 1
                           r2 = dist2_pbc(x_crs(:, p1), x_crs(:, p2), box_len, box_half)
                           if (r2 <= rp2) call dsu_union_serial(parent_size, p1, p2)
                        end do
                     end do
                     ! 收集去重 jcell > icell
                     num_nbr = 0
                     do dz = -1, 1
                        if (.not. periodic_dim(3) .and. (iz + dz < 1 .or. iz + dz > Nz)) cycle
                        jz = modulo(iz + dz - 1, Nz) + 1
                        do dy = -1, 1
                           if (.not. periodic_dim(2) .and. (iy + dy < 1 .or. iy + dy > Ny)) cycle
                           jy = modulo(iy + dy - 1, Ny) + 1
                           do dx_step = -1, 1
                              if (.not. periodic_dim(1) .and. (ix + dx_step < 1 .or. ix + dx_step > Nx)) cycle
                              jx = modulo(ix + dx_step - 1, Nx) + 1
                              jcell = jx + Nx * (jy - 1) + Nx * Ny * (jz - 1)
                              if (jcell <= icell) cycle
                              duplicate = .false.
                              do k = 1, num_nbr
                                 if (nbr_cells(k) == jcell) then
                                    duplicate = .true.
                                    exit
                                 end if
                              end do
                              if (.not. duplicate) then
                                 num_nbr = num_nbr + 1
                                 nbr_cells(num_nbr) = jcell
                              end if
                           end do
                        end do
                     end do
                     do k = 1, num_nbr
                        jcell = nbr_cells(k)
                        if (cell_ptr(jcell) >= cell_ptr(jcell + 1)) cycle
                        do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                           do p2 = cell_ptr(jcell), cell_ptr(jcell + 1) - 1
                              r2 = dist2_pbc(x_crs(:, p1), x_crs(:, p2), box_len, box_half)
                              if (r2 <= rp2) call dsu_union_serial(parent_size, p1, p2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do

         else
            allocate(layer_done(Nz), stat=ios)
            if (ios /= 0) then
               alloc_err = 3
               deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
               deallocate(x_crs, v_crs, mp_crs)
               return
            end if
            layer_done = .false.
            !$omp parallel num_threads(nest_threads) default(shared) &
            !$omp private(ithread, my_start, step, iz, claimed, iy, ix, icell, p1, p2, r2, &
            !$omp         dz, dy, dx_step, jx, jy, jz, jcell, num_nbr, nbr_cells, duplicate, k)
            ithread = omp_get_thread_num() + 1
            my_start = int(int(Nz, 8) * (ithread - 1) / omp_get_num_threads()) + 1
            do step = 0, Nz - 1
               iz = mod(my_start - 1 + step, Nz) + 1
               !$omp atomic capture
               claimed = layer_done(iz)
               layer_done(iz) = .true.
               !$omp end atomic
               if (claimed) cycle

               do iy = 1, Ny
                  do ix = 1, Nx
                     icell = ix + Nx * (iy - 1) + Nx * Ny * (iz - 1)
                     if (cell_ptr(icell) >= cell_ptr(icell + 1)) cycle
                     ! 本 cell 内上三角
                     do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                        do p2 = p1 + 1, cell_ptr(icell + 1) - 1
                           r2 = dist2_pbc(x_crs(:, p1), x_crs(:, p2), box_len, box_half)
                           if (r2 <= rp2) call fof_atomic_union(parent_size, p1, p2)
                        end do
                     end do
                     ! 收集去重 jcell > icell
                     num_nbr = 0
                     do dz = -1, 1
                        if (.not. periodic_dim(3) .and. (iz + dz < 1 .or. iz + dz > Nz)) cycle
                        jz = modulo(iz + dz - 1, Nz) + 1
                        do dy = -1, 1
                           if (.not. periodic_dim(2) .and. (iy + dy < 1 .or. iy + dy > Ny)) cycle
                           jy = modulo(iy + dy - 1, Ny) + 1
                           do dx_step = -1, 1
                              if (.not. periodic_dim(1) .and. (ix + dx_step < 1 .or. ix + dx_step > Nx)) cycle
                              jx = modulo(ix + dx_step - 1, Nx) + 1
                              jcell = jx + Nx * (jy - 1) + Nx * Ny * (jz - 1)
                              if (jcell <= icell) cycle
                              duplicate = .false.
                              do k = 1, num_nbr
                                 if (nbr_cells(k) == jcell) then
                                    duplicate = .true.
                                    exit
                                 end if
                              end do
                              if (.not. duplicate) then
                                 num_nbr = num_nbr + 1
                                 nbr_cells(num_nbr) = jcell
                              end if
                           end do
                        end do
                     end do
                     do k = 1, num_nbr
                        jcell = nbr_cells(k)
                        if (cell_ptr(jcell) >= cell_ptr(jcell + 1)) cycle
                        do p1 = cell_ptr(icell), cell_ptr(icell + 1) - 1
                           do p2 = cell_ptr(jcell), cell_ptr(jcell + 1) - 1
                              r2 = dist2_pbc(x_crs(:, p1), x_crs(:, p2), box_len, box_half)
                              if (r2 <= rp2) call fof_atomic_union(parent_size, p1, p2)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
            !$omp end parallel
            deallocate(layer_done)
         end if
      end if

      ! 统计各分量粒子数并筛选候选
      allocate(root_size(np_in), root_map(np_in), stat=ios)
      if (ios /= 0) then
         alloc_err = 3
         deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
         deallocate(x_crs, v_crs, mp_crs)
         return
      end if
      root_size = 0
      root_map = 0

      do p = 1, np_in
         root_p = dsu_find_root(parent_size, p)
         root_size(root_p) = root_size(root_p) + 1
      end do

      cand_cnt = 0
      mem_cnt = 0
      do root_p = 1, np_in
         if (root_size(root_p) >= np_min) then
            cand_cnt = cand_cnt + 1
            root_map(root_p) = cand_cnt
            mem_cnt = mem_cnt + root_size(root_p)
         end if
      end do

      cand_out%ncand = cand_cnt
      if (cand_cnt > 0) then
         allocate(cand_out%cand_np(cand_cnt), cand_out%cand_head(cand_cnt), stat=ios)
         if (ios == 0) allocate(cand_out%cand_members(mem_cnt), cand_out%member_ptr(cand_cnt + 1), stat=ios)
         if (ios == 0) allocate(cand_out%xm(3, cand_cnt), cand_out%vm(3, cand_cnt), stat=ios)
         if (ios /= 0) then
            alloc_err = 3
            if (allocated(cand_out%cand_np)) deallocate(cand_out%cand_np, cand_out%cand_head)
            if (allocated(cand_out%cand_members)) deallocate(cand_out%cand_members, cand_out%member_ptr)
            if (allocated(cand_out%xm)) deallocate(cand_out%xm, cand_out%vm)
            cand_out%ncand = 0
            deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
            deallocate(x_crs, v_crs, mp_crs, root_size, root_map)
            return
         end if

         cand_cnt = 0
         cand_out%member_ptr(1) = 1
         do root_p = 1, np_in
            if (root_size(root_p) >= np_min) then
               cand_cnt = cand_cnt + 1
               cand_out%cand_np(cand_cnt) = root_size(root_p)
               cand_out%member_ptr(cand_cnt + 1) = cand_out%member_ptr(cand_cnt) + root_size(root_p)
            end if
         end do

         ! 将成员写入 cand_members
         allocate(cand_cursor(cand_cnt), stat=ios)
         if (ios /= 0) then
            alloc_err = 3
            deallocate(cand_out%cand_np, cand_out%cand_head, cand_out%cand_members, cand_out%member_ptr, cand_out%xm, cand_out%vm)
            cand_out%ncand = 0
            deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
            deallocate(x_crs, v_crs, mp_crs, root_size, root_map)
            return
         end if
         cand_cursor = cand_out%member_ptr(1:cand_cnt)

         do p = 1, np_in
            root_p = dsu_find_root(parent_size, p)
            g = root_map(root_p)
            if (g > 0) then
               cand_out%cand_members(cand_cursor(g)) = pid_reg(p)
               cand_cursor(g) = cand_cursor(g) + 1
            end if
         end do
         deallocate(cand_cursor)

         ! 计算每个候选的初步质心与均速（以团块首成员为锚点做周期展开）
         do g = 1, cand_cnt
            anchor_idx = cand_out%cand_members(cand_out%member_ptr(g))
            xp_anchor_cand = xp_halo(:, anchor_idx)
            w_sum = 0.0d0
            qc_d = 0.0d0
            vc_d = 0.0d0
            do p = cand_out%member_ptr(g), cand_out%member_ptr(g + 1) - 1
               orig_idx = cand_out%cand_members(p)
               wi = mp_halo(orig_idx)
               w_sum = w_sum + wi
               dx = modulo(xp_halo(:, orig_idx) - xp_anchor_cand + box_half, box_len) - box_half
               qc_d = qc_d + wi * dx
               vc_d = vc_d + wi * vp_halo(:, orig_idx)
            end do
            if (w_sum > 0.0d0) then
               cand_out%xm(:, g) = real(qc_d / w_sum, 4)
               cand_out%vm(:, g) = real(vc_d / w_sum, 4)
            else
               cand_out%xm(:, g) = 0.0
               cand_out%vm(:, g) = 0.0
            end if
         end do
      end if

      deallocate(cell_ptr, cell_cursor, cell_cursor_save, parent_size, pid_reg)
      deallocate(x_crs, v_crs, mp_crs, root_size, root_map)
   end subroutine process_single_region

end subroutine fof_fortran
