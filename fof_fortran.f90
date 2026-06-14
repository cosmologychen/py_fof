!===============================================================================
! CUBE_FoF: Friends-of-Friends (FOF) halo finder for cosmological simulations
!
! 功能说明:
!   对给定的粒子坐标执行FOF算法，识别满足链接长度条件的粒子团块
!   支持周期性边界条件、束缚检查（基于动能与势能）
!
! 作者: [Your Name]
! 日期: 2024
!
! 输入参数:
!   xp_halo      : 粒子坐标数组，shape=(3, total_np)，单位：格点
!   vp_halo      : 粒子速度数组，shape=(3, total_np)，单位：格点/时间步
!   np_halo      : 每个halo的粒子数量数组，shape=(nhalo,)
!   nhalo        : halo数量
!   L_Box        : 盒子边长，用于周期性边界条件，单位：格点
!   linking_length: FOF链接长度，单位：格点
!   np_halo_min  : 最小halo粒子数阈值
!   bound_check  : 是否执行束缚检查
!
! 输出参数:
!   in_halo      : 逻辑数组，标记每个粒子是否属于有效halo
!
! 编译方法:
!   python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags="-O3"
!===============================================================================

subroutine fof_fortran(xp_halo, vp_halo, np_halo, nhalo, L_Box, &
   linking_length, np_halo_min, bound_check, verbose, in_halo, total_np, ang_mom)
   implicit none

   !---------------------------------------------------------------------------
   ! 输入参数声明
   !---------------------------------------------------------------------------
   integer, intent(in) :: nhalo                          ! halo数量
   integer, intent(in) :: total_np                       ! 总粒子数
   integer, intent(in) :: np_halo(nhalo)                 ! 每个halo的粒子数
   real(8), intent(in) :: xp_halo(3, total_np)           ! 粒子坐标
   real(8), intent(in) :: vp_halo(3, total_np)           ! 粒子速度
   real(8), intent(in) :: L_Box                          ! 盒子边长
   real(8), intent(in) :: linking_length                 ! 链接长度
   integer, intent(in) :: np_halo_min                    ! 最小粒子数阈值
   integer, intent(in) :: bound_check                    ! 是否束缚检查 (0=False, 1=True)
   integer, intent(in) :: verbose                        ! 详细输出级别 (0=关闭, 1=基本, 2=详细)
   integer, intent(out) :: in_halo(total_np)             ! 输出：粒子是否在halo中 (0/1)
   real(8), intent(out) :: ang_mom(3, nhalo)             ! 输出：每个halo的角动量

   !f2py intent(in) :: nhalo
   !f2py intent(in) :: total_np
   !f2py intent(in) :: np_halo
   !f2py intent(in) :: xp_halo
   !f2py intent(in) :: vp_halo
   !f2py intent(in) :: L_Box
   !f2py intent(in) :: linking_length
   !f2py intent(in) :: np_halo_min
   !f2py intent(in) :: bound_check
   !f2py intent(in) :: verbose
   !f2py intent(out) :: in_halo
   !f2py intent(out) :: ang_mom
   !f2py depend(nhalo) :: np_halo, ang_mom
   !f2py depend(total_np) :: xp_halo, vp_halo, in_halo

   !---------------------------------------------------------------------------
   ! 本地变量声明
   !---------------------------------------------------------------------------
   real(8), allocatable :: xp_local(:,:)    ! 本地坐标副本

   ! 链表相关变量
   integer, allocatable :: llgp(:)          ! 粒子链表
   integer, allocatable :: hcgp(:)          ! 链头
   integer, allocatable :: ecgp(:)          ! 链尾

   ! 网格链表变量
   integer, allocatable :: hoc_local(:,:,:) ! 网格头粒子
   integer, allocatable :: ll_local(:)      ! 局部链表

   ! 索引数组
   integer, allocatable :: i_start_array(:) ! 每个halo的起始索引

   ! 坐标相关变量
   real(8) :: xp_range(3)
   real(8) :: xp_min(3)
   real(8) :: xp_mean(3)
   real(8) :: xp_zero(3)
   real(8) :: L_mesh
   real(8) :: dxv(6)

   ! 循环和计数变量
   integer :: np_max
   integer :: ihalo
   integer :: np
   integer :: istart
   integer :: ip, jp, jp_prev
   integer :: N_x, N_y, N_z
   integer :: iq1, iq2, iq3
   integer :: i
   integer :: np_fof
   integer :: np_bound
   integer :: iteration

   ! 距离和速度计算
   real(8) :: rsq
   real(8) :: dx(3)
   real(8) :: dv(3)
   real(8) :: dx_rel(3)
   real(8) :: dv_rel(3)
   real(8) :: r
   real(8) :: v_rel_sq
   real(8) :: L_Box2
   real(8) :: rp2

   ! 网格索引
   integer :: idx3(3)
   integer :: jq(3)
   integer :: i_neighbor

   ! 邻居偏移数组
   integer, parameter :: N_neighbors = 13
   integer :: ijk(3, N_neighbors)

   !---------------------------------------------------------------------------
   ! 初始化邻居偏移数组
   !---------------------------------------------------------------------------
   ijk(:, 1)  = [-1, -1, -1];  ijk(:, 2)  = [-1, -1,  0]
   ijk(:, 3)  = [-1, -1,  1];  ijk(:, 4)  = [-1,  0, -1]
   ijk(:, 5)  = [-1,  0,  0];  ijk(:, 6)  = [-1,  0,  1]
   ijk(:, 7)  = [-1,  1, -1];  ijk(:, 8)  = [-1,  1,  0]
   ijk(:, 9)  = [-1,  1,  1];  ijk(:, 10) = [ 0, -1, -1]
   ijk(:, 11) = [ 0, -1,  0];  ijk(:, 12) = [ 0, -1,  1]
   ijk(:, 13) = [ 0,  0, -1]

   !---------------------------------------------------------------------------
   ! 初始化常量
   !---------------------------------------------------------------------------
   L_Box2 = L_Box / 2.0d0
   rp2 = linking_length * linking_length

   !---------------------------------------------------------------------------
   ! 创建坐标的本地副本（因为需要修改坐标）
   !---------------------------------------------------------------------------
   np_max = total_np
   allocate(xp_local(3, np_max))
   xp_local = xp_halo

   !---------------------------------------------------------------------------
   ! 计算每个halo的起始索引
   !---------------------------------------------------------------------------
   allocate(i_start_array(nhalo))
   i_start_array(1) = 1
   do ihalo = 2, nhalo
      i_start_array(ihalo) = i_start_array(ihalo-1) + np_halo(ihalo-1)
   enddo

   !---------------------------------------------------------------------------
   ! 初始化全局链表结构
   !---------------------------------------------------------------------------
   allocate(llgp(np_max), hcgp(np_max), ecgp(np_max))

   do ip = 1, np_max
      hcgp(ip) = ip
   enddo
   llgp = 0
   ecgp = 0
   in_halo = 0
   ang_mom = 0.0d0

   !---------------------------------------------------------------------------
   ! 主循环：遍历每个halo
   !---------------------------------------------------------------------------
   do ihalo = 1, nhalo
      istart = i_start_array(ihalo)
      np = np_halo(ihalo)

      ! 输出halo处理进度
      if (verbose >= 1) then
         write(*, '(A,I4,A,I6,A)') '  Processing halo ', ihalo, '/', nhalo
         write(*, '(A,I6,A)') '    N_particles: ', np
      endif

      !-----------------------------------------------------------------------
      ! 计算质心坐标
      !-----------------------------------------------------------------------
      dxv(1:3) = 0.0d0
      do ip = istart + 1, istart + np - 1
         dxv(1:3) = dxv(1:3) + modulo((xp_local(:, ip) - xp_local(:, istart)) + L_Box2, L_Box) - L_Box2
      enddo
      xp_mean = modulo(dxv(1:3) / dble(np) + xp_local(:, istart), L_Box)
      do ip = istart, istart + np - 1
         xp_local(:, ip) = xp_local(:, ip) - xp_mean
      enddo

      !-----------------------------------------------------------------------
      ! 构建局部网格链表
      !-----------------------------------------------------------------------
      xp_min = minval(xp_local(:, istart:istart+np-1), dim=2)
      xp_range = maxval(xp_local(:, istart:istart+np-1), dim=2) - xp_min

      L_mesh = max(1.0d0, maxval(xp_range)**3 / dble(np))

      N_x = max(1, ceiling(xp_range(1) / L_mesh))
      N_y = max(1, ceiling(xp_range(2) / L_mesh))
      N_z = max(1, ceiling(xp_range(3) / L_mesh))

      ! 输出网格大小
      if (verbose >= 2) then
         write(*, '(A,I4,A,I4,A,I4)') '    N_xyz: ', N_x, ' x', N_y, ' x', N_z
      endif


      allocate(hoc_local(N_x, N_y, N_z), ll_local(np))
      hoc_local = 0
      ll_local = 0

      do i = 1, np
         idx3 = floor((xp_local(:, istart+i-1) - xp_min) / L_mesh) + 1
         idx3(1) = max(1, min(idx3(1), N_x))
         idx3(2) = max(1, min(idx3(2), N_y))
         idx3(3) = max(1, min(idx3(3), N_z))
         ll_local(i) = hoc_local(idx3(1), idx3(2), idx3(3))
         hoc_local(idx3(1), idx3(2), idx3(3)) = i
      enddo

      !-----------------------------------------------------------------------
      ! FOF搜索
      !-----------------------------------------------------------------------
      do iq3 = 1, N_z
         do iq2 = 1, N_y
            do iq1 = 1, N_x
               ip = hoc_local(iq1, iq2, iq3)

               do while (ip /= 0)
                  jp = ll_local(ip)
                  do while (jp /= 0)
                     rsq = sum((xp_local(1:3, ip+istart-1) - xp_local(1:3, jp+istart-1))**2)
                     if (rsq <= rp2) call merge_chain(ip+istart-1, jp+istart-1)
                     jp = ll_local(jp)
                  enddo

                  do i_neighbor = 1, N_neighbors
                     jq(1) = modulo(iq1 + ijk(1, i_neighbor)-1, N_x)+1
                     jq(2) = modulo(iq2 + ijk(2, i_neighbor)-1, N_y)+1
                     jq(3) = modulo(iq3 + ijk(3, i_neighbor)-1, N_z)+1

                     jp = hoc_local(jq(1), jq(2), jq(3))
                     do while (jp /= 0)
                        rsq = sum((xp_local(1:3, ip+istart-1) - xp_local(1:3, jp+istart-1))**2)
                        if (rsq <= rp2) call merge_chain(ip+istart-1, jp+istart-1)
                        jp = ll_local(jp)
                     enddo
                  enddo

                  ip = ll_local(ip)
               enddo
            enddo
         enddo
      enddo

      deallocate(hoc_local, ll_local)

      !-----------------------------------------------------------------------
      ! 统计FOF链并标记有效halo
      !-----------------------------------------------------------------------
      do i = istart, istart + np - 1
         if (hcgp(i) /= i .and. hcgp(i) /= 0) then
            ip = hcgp(i)

            np_fof = 0
            jp = ip
            do while (jp /= 0)
               np_fof = np_fof + 1
               jp = llgp(jp)
            enddo

            if (np_fof >= np_halo_min) then
               dxv(1:3) = 0.0d0
               jp = ip
               do while (jp /= 0)
                  dxv(1:3) = dxv(1:3) + xp_local(:, jp)
                  jp = llgp(jp)
               enddo
               dxv(1:3) = modulo(dxv(1:3) / dble(np_fof) + xp_mean, L_Box)

               if (bound_check == 1) then
                  iteration = 0
                  do while (iteration < 3)
                     iteration = iteration + 1
                     np_bound = 0
                     jp_prev = 0
                     jp = ip

                     ! 输出bound_check轮次
                     if (verbose >= 1) then
                        write(*, '(A,I1,A,I6)') '      Bound check iteration ', iteration, &
                           ', N_particles before: ', np_fof
                     endif

                     dxv(4:6) = 0.0d0
                     jp = ip
                     do while (jp /= 0)
                        dxv(4:6) = dxv(4:6) + vp_halo(:, jp)
                        jp = llgp(jp)
                     enddo
                     dxv(4:6) = dxv(4:6) / dble(np_fof)

                     jp = ip
                     do while (jp /= 0)
                        dx = modulo((xp_local(1:3, jp) - dxv(1:3)) + L_Box2, L_Box) - L_Box2
                        dv = vp_halo(1:3, jp) - dxv(4:6)
                        r = norm2(dx)
                        v_rel_sq = sum(dv**2)

                        if (r > 0.0d0 .and. 0.5d0 * v_rel_sq < dble(np_fof) / r) then
                           np_bound = np_bound + 1
                           jp_prev = jp
                           jp = llgp(jp)
                        else
                           if (jp == ip) then
                              ip = llgp(jp)
                              hcgp(i) = ip
                              if (ip /= 0) hcgp(ip) = ip
                              jp = ip
                           else
                              llgp(jp_prev) = llgp(jp)
                              jp = llgp(jp_prev)
                           endif
                        endif
                     enddo

                     ! 输出丢弃了多少粒子
                     if (verbose >= 1) then
                        write(*, '(A,I6,A,I6)') '      N_particles after: ', np_bound, &
                           ', removed: ', np_fof - np_bound
                     endif

                     if (np_bound == np_fof) exit
                     if (np_bound < np_halo_min) exit

                     np_fof = np_bound
                     dxv(1:3) = 0.0d0
                     jp = ip
                     do while (jp /= 0)
                        dxv(1:3) = dxv(1:3) + xp_local(:, jp)
                        jp = llgp(jp)
                     enddo
                     dxv(1:3) = modulo(dxv(1:3) / dble(np_fof) + xp_mean, L_Box)
                  enddo
               endif

               if (np_fof >= np_halo_min) then
                  if (bound_check /= 1) then
                     dxv(4:6) = 0.0d0
                     jp = ip
                     do while (jp /= 0)
                        dxv(4:6) = dxv(4:6) + vp_halo(:, jp)
                        jp = llgp(jp)
                     enddo
                     dxv(4:6) = dxv(4:6) / dble(np_fof)
                  endif

                  jp = ip
                  do while (jp /= 0)
                     in_halo(jp) = 1
                     dx_rel = modulo((xp_local(1:3, jp) - dxv(1:3)) + L_Box2, L_Box) - L_Box2
                     dv_rel = vp_halo(1:3, jp) - dxv(4:6)
                     ang_mom(1, ihalo) = ang_mom(1, ihalo) + (dx_rel(2)*dv_rel(3) - dx_rel(3)*dv_rel(2))
                     ang_mom(2, ihalo) = ang_mom(2, ihalo) + (dx_rel(3)*dv_rel(1) - dx_rel(1)*dv_rel(3))
                     ang_mom(3, ihalo) = ang_mom(3, ihalo) + (dx_rel(1)*dv_rel(2) - dx_rel(2)*dv_rel(1))
                     jp = llgp(jp)
                  enddo
               endif
            endif
         endif
      enddo
   enddo




   !---------------------------------------------------------------------------
   ! 释放内存
   !---------------------------------------------------------------------------
   deallocate(xp_local, llgp, hcgp, ecgp, i_start_array)

contains

   !---------------------------------------------------------------------------
   ! merge_chain: 合并两个FOF链
   !---------------------------------------------------------------------------
   subroutine merge_chain(ii, jj)
      integer, intent(in) :: ii, jj
      integer :: ihead, jhead, iend, jend, ipart

      jend = merge(jj, ecgp(jj), ecgp(jj) == 0)
      iend = merge(ii, ecgp(ii), ecgp(ii) == 0)

      if (iend == jend) return

      ihead = max(hcgp(ii), hcgp(iend))
      jhead = max(hcgp(jj), hcgp(jend))

      ipart = jhead
      do while (ipart /= 0)
         ecgp(ipart) = iend
         ipart = llgp(ipart)
      enddo

      llgp(jend) = ihead
      ecgp(ii) = iend
      hcgp(iend) = jhead
      hcgp(jend) = 0
   end subroutine merge_chain

end subroutine fof_fortran
