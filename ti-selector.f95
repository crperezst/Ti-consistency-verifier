	
	program main
	
	character :: comp_file*40,image_file*40,output_file*40,parameters*60,par*60,empty*40
	
	integer :: comp_num_in_file,blocks_x_ti,blocks_y_ti,blocks_z_ti
	integer :: search_rad_x,search_rad_y,search_rad_z
	integer :: min_conditioning_nodes,max_conditioning_nodes,random_seed
	integer :: blocks_in_ti,integer_random,nvar,ti_num,ti_facie
	integer :: comp_out_of_bounds,migrated_nodes,repeated_nodes
	integer :: xk,yk,zk
	integer :: blocks_in_mti,blocks_in_mcomp_grid
	integer :: blocks_x_mti,blocks_y_mti,blocks_z_mti
	integer :: blocks_x_mcomp_grid,blocks_y_mcomp_grid,blocks_z_mcomp_grid
	integer :: xt,yt,zt,nt,mid_n_tmp
	integer :: rinteger
	integer :: offx,offy,offz
	integer :: startx,starty,startz
	integer :: coordtempx,coordtempy,coordtempz
	integer :: visited_node,cond_node
	integer :: valid_events,invalid_events,total_events
	integer :: min_equal_nodes
	integer :: tioffx,tioffy,tioffz,tiofffacie
	integer :: coordtix,coordtiy,coordtiz
	integer :: equal_nodes
	integer :: rintegerx,rintegery,rintegerz
	integer :: blocks_x_comp,blocks_y_comp,blocks_z_comp
	integer :: blocks_in_comp,blocks_to_check
	
	integer :: aux1,aux3,aux2
	integer :: seed_counter,read_file_counter
	integer :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
	integer :: c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28
	integer :: c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,c41
	integer :: c42,c43,c44,c45,c46,c47,c48,c49,c50,c51,c52
	integer :: checker
	integer :: mode
	integer :: column
	integer :: var_type,acc_type
	
	real :: auxr1,auxr2
	
	real :: ti_x_block_size,ti_y_block_size,ti_z_block_size
	real :: image_min_x_centroid,image_min_y_centroid,image_min_z_centroid
	real :: cond_event_tolerance,random,randomx,randomy,randomz
	real :: image_min_x,image_min_y,image_min_z
	real :: image_max_x,image_max_y,image_max_z
	real :: xcomp,ycomp,zcomp,fcomp
	real :: xi,yi,zi,xj,yj,zj
	real :: orderer
	real :: comp_x_block_size,comp_y_block_size,block_z_comp_size
	real :: comp_min_x_centroid,comp_min_y_centroid,comp_min_z_centroid
	real :: ti_check_proportion
	real :: comp_min_x,comp_min_y,comp_min_z,comp_max_x,comp_max_y,comp_max_z
	real :: acc_thresh
	
	
	real,allocatable :: comp(:,:)
	real,allocatable :: ti(:,:,:,:)
	real,allocatable :: tif(:)
	real,allocatable :: mti(:,:,:,:)
	real,allocatable :: comp_grid(:,:,:)
	real,allocatable :: mcomp_grid(:,:,:)
	real,allocatable :: tmp(:,:)
	real,allocatable :: tmp_ord(:,:)
	real,allocatable :: c_event(:,:)
	real,allocatable :: cm(:,:)
	real,allocatable :: ncm(:,:)
	real,allocatable :: cm2(:,:)
	logical :: testfl
	real :: absolutecomp
	
	par='ti-selector.par'
	empty=''
	
	write(*,*) 'TI-Selector 1.0'
	write(*,*) 'Which parameter file do you want to use?'
	read (*,'(a40)') parameters
	if (parameters.eq.empty) then
		inquire(file=par,exist=testfl)
		if(.not.testfl) then
			write(*,*) 'ERROR - the parameter file does not exist,'
			write(*,*) '        writing a new parameters file  '
			open(500,file=par)
			write(500,*) '                                   Parameters'
			write(500,*) '                                  ************'
			write(500,*) ''
			write(500,*) 'START OF PARAMETERS:'
			write(500,*) 'samples.out                                       - samples datafile name'
			write(500,*) '1091                                              - sample number'
			write(500,*) '100 0.5 1                                         - nx, xmn, xsiz; in sample migration grid'
			write(500,*) '100 0.5 1                                         - ny, ymn, ysiz; in sample migration grid'
			write(500,*) '1   0.5 1                                         - nz, zmn, zsiz; in sample migration grid'
			write(500,*) 'ti-set.out                                        - training images datafile name'
			write(500,*) '200 0.5 1                                         - nx, xmn, xsiz; in training images'
			write(500,*) '200 0.5 1                                         - ny, ymn, ysiz; in training images'
			write(500,*) '1   0.5 1                                         - nz, zmn, zsiz; in training images'
			write(500,*) '25 25 0                                           - search radii in x,y and z directions'
			write(500,*) '15                                                - minimum event order to analyze'
			write(500,*) '25                                                - maximum event order to analyze'
			write(500,*) '0                                                 - compatibility mode: 0 for relative, 1 for absolute'
			write(500,*) '1                                                 - column (for absolute mode)'
			write(500,*) '0                                                 - tolerance parameter'
			write(500,*) '0                                                 - variable type'
			write(500,*) '0.1                                               - acceptance threshold (continuous variables)'
			write(500,*) '0                                                 - training image fraction to check'
			write(500,*) '12345                                             - random seed'
			write(500,*) 'output.out                                        - output file name'
			close(500)
		stop
		end if
	else
		inquire(file=parameters,exist=testfl)
		if(.not.testfl) then
			write(*,*) 'ERROR - the parameter file does not exist,'
			write(*,*) '        check for the file and try again'
			stop
		end if
	end if
	
	if (parameters.eq.empty) then
		open(20,file=par)
		
		read(20,*,err=98)
		read(20,*,err=98)
		read(20,*,err=98)
		read(20,*,err=98)
			
		read(20,'(a40)',err=98) comp_file
		read(20,*,err=98) comp_num_in_file
		read(20,*,err=98) blocks_x_comp,comp_min_x_centroid,comp_x_block_size
		read(20,*,err=98) blocks_y_comp,comp_min_y_centroid,comp_y_block_size
		read(20,*,err=98) blocks_z_comp,comp_min_z_centroid,comp_z_block_size
		read(20,'(a40)',err=98) image_file
		read(20,*,err=98) blocks_x_ti,image_min_x_centroid,ti_x_block_size
		read(20,*,err=98) blocks_y_ti,image_min_y_centroid,ti_y_block_size
		read(20,*,err=98) blocks_z_ti,image_min_z_centroid,ti_z_block_size
		read(20,*,err=98) search_rad_x,search_rad_y,search_rad_z
		read(20,*,err=98) min_conditioning_nodes
		read(20,*,err=98) max_conditioning_nodes
		read(20,*,err=98) mode
		read(20,*,err=98) column
		read(20,*,err=98) cond_event_tolerance
		read(20,*,err=98) var_type
		read(20,*,err=98) acc_thresh
		read(20,*,err=98) ti_check_proportion
		read(20,*,err=98) random_seed
		read(20,'(a40)',err=98) output_file
	else
		open(10,file=parameters)
	
		read(10,*,err=98)
		read(10,*,err=98)
		read(10,*,err=98)
		read(10,*,err=98)
	
		read(10,'(a40)',err=98) comp_file
		read(10,*,err=98) comp_num_in_file
		read(10,*,err=98) blocks_x_comp,comp_min_x_centroid,comp_x_block_size
		read(10,*,err=98) blocks_y_comp,comp_min_y_centroid,comp_y_block_size
		read(10,*,err=98) blocks_z_comp,comp_min_z_centroid,comp_z_block_size
		read(10,'(a40)',err=98) image_file
		read(10,*,err=98) blocks_x_ti,image_min_x_centroid,ti_x_block_size
		read(10,*,err=98) blocks_y_ti,image_min_y_centroid,ti_y_block_size
		read(10,*,err=98) blocks_z_ti,image_min_z_centroid,ti_z_block_size
		read(10,*,err=98) search_rad_x,search_rad_y,search_rad_z
		read(10,*,err=98) min_conditioning_nodes
		read(10,*,err=98) max_conditioning_nodes
		read(10,*,err=98) mode
		read(10,*,err=98) column
		read(10,*,err=98) cond_event_tolerance
		read(10,*,err=98) var_type
		read(10,*,err=98) acc_thresh
		read(10,*,err=98) ti_check_proportion
		read(10,*,err=98) random_seed
		read(10,'(a40)',err=98) output_file
	end if
	
	allocate (comp(comp_num_in_file,4))
	allocate (c_event(max_conditioning_nodes,4))
	
	if (mode.ne.1.and.mode.ne.0) then
		goto 95
	end if
	
	comp_file=adjustl(comp_file)
	image_file=adjustl(image_file)
	output_file=adjustl(output_file)
	
	inquire(file=comp_file,exist=testfl)
	if(.not.testfl) then
		write(*,*) 'ERROR - the composites file does not exist,'
		write(*,*) '        check for the file and try again  '
		stop
	endif
	
	inquire(file=image_file,exist=testfl)
	if(.not.testfl) then
		write(*,*) 'ERROR - the trainig images file does not exist,'
		write(*,*) '        check for the file and try again  '
		stop
	endif
	
	if (mode==1) then
		ti_check_proportion=1
	end if
	
	if (mode==0.and.ti_check_proportion==0) then
		comp_out_of_bounds=0
		migrated_nodes=0
		repeated_nodes=0
		aux4=0
		
		
		blocks_in_ti=blocks_x_ti*blocks_y_ti*blocks_z_ti
		blocks_in_comp=blocks_x_comp*blocks_y_comp*blocks_z_comp
		image_min_x=image_min_x_centroid-(ti_x_block_size/2)
		image_min_y=image_min_y_centroid-(ti_y_block_size/2)
		image_min_z=image_min_z_centroid-(ti_z_block_size/2)
		comp_min_x=comp_min_x_centroid-(comp_x_block_size/2)
		comp_min_y=comp_min_y_centroid-(comp_y_block_size/2)
		comp_min_z=comp_min_z_centroid-(comp_z_block_size/2)
		image_max_x=image_min_x_centroid+(ti_x_block_size*(blocks_x_ti-1))+(ti_x_block_size/2)
		image_max_y=image_min_y_centroid+(ti_y_block_size*(blocks_y_ti-1))+(ti_y_block_size/2)
		image_max_z=image_min_z_centroid+(ti_z_block_size*(blocks_z_ti-1))+(ti_z_block_size/2)
		comp_max_x=comp_min_x_centroid+(comp_x_block_size*(blocks_x_comp-1))+(comp_x_block_size/2)
		comp_max_y=comp_min_y_centroid+(comp_y_block_size*(blocks_y_comp-1))+(comp_y_block_size/2)
		comp_max_z=comp_min_z_centroid+(comp_z_block_size*(blocks_z_comp-1))+(comp_z_block_size/2)
		blocks_x_mti=blocks_x_ti+2*search_rad_x
		blocks_y_mti=blocks_y_ti+2*search_rad_y
		blocks_z_mti=blocks_z_ti+2*search_rad_z
		blocks_in_mti=blocks_x_mti*blocks_y_mti*blocks_z_mti
		blocks_x_mcomp_grid=blocks_x_comp+2*search_rad_x
		blocks_y_mcomp_grid=blocks_y_comp+2*search_rad_y
		blocks_z_mcomp_grid=blocks_z_comp+2*search_rad_z
		blocks_in_mcomp_grid=blocks_x_mcomp_grid*blocks_y_mcomp_grid*blocks_z_mcomp_grid
		xt=search_rad_x*2+1
		yt=search_rad_y*2+1
		zt=search_rad_z*2+1
		nt=xt*yt*zt
		auxr1=nt/2
		mid_n_tmp=floor(auxr1)+1
		
		
		call init_random_seed()
		
		do seed_counter=1,random_seed
			call random_number(random)
		end do
		
		open(1,file=comp_file)
		read(1,*,err=99)
		read(1,*,err=99) nvar
		do read_file_counter=1,nvar
			read(1,*,err=99)
		end do
		
		
		open(2,file=image_file)
		read(2,*,err=97)
		read(2,*,err=97) ti_num
		do read_file_counter=1,ti_num
			read(2,*,err=97)
		end do

		allocate (ti(blocks_x_ti,blocks_y_ti,blocks_z_ti,ti_num))
		allocate (tif(ti_num))
		allocate (mti(blocks_x_mti,blocks_y_mti,blocks_z_mti,ti_num))
		allocate (comp_grid(blocks_x_comp,blocks_y_comp,blocks_z_comp))
		allocate (mcomp_grid(blocks_x_mcomp_grid,blocks_y_mcomp_grid,blocks_z_mcomp_grid))
		allocate (tmp(nt,5))
		allocate (tmp_ord(nt,4))
		allocate (cm(1,ti_num))
		allocate (ncm(1,ti_num))
		allocate (cm2(1,ti_num))
		
		do c40=1,ti_num
			cm(1,c40)=0
			ncm(1,c40)=0
			cm2(1,c40)=0
		end do
		
		open(3,file=output_file)
		
		do c1=1,comp_num_in_file
			read(1,*,err=99) xcomp,ycomp,zcomp,fcomp
			comp(c1,1)=xcomp
			comp(c1,2)=ycomp
			comp(c1,3)=zcomp
			comp(c1,4)=fcomp
		end do
		
		do c2=1,blocks_z_ti
			do c3=1,blocks_y_ti
				do c4=1,blocks_x_ti
					read(2,*,err=97) (tif(c49),c49=1,ti_num)
					do c50=1,ti_num
						ti(c4,c3,c2,c50)=tif(c50)
					end do
				end do
			end do
		end do
		
		do c5=1,blocks_z_comp
			do c6=1,blocks_y_comp
				do c7=1,blocks_x_comp
					comp_grid(c7,c6,c5)=-111
				end do
			end do
		end do
		
		do c8=1,comp_num_in_file
			if (comp(c8,1)<comp_min_x .or. comp(c8,1)>comp_max_x .or. comp(c8,2)<comp_min_y .or. comp(c8,2)>comp_max_y) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			if (comp(c8,3)<comp_min_z .or. comp(c8,3)>comp_max_z) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			xi=comp(c8,1)-comp_min_x
			yi=comp(c8,2)-comp_min_y
			zi=comp(c8,3)-comp_min_z
			xj=xi/comp_x_block_size+1
			yj=yi/comp_y_block_size+1
			zj=zi/comp_z_block_size+1
			xk=floor(xj)
			yk=floor(yj)
			zk=floor(zj)
			if (comp_grid(xk,yk,zk)==-111) then
				migrated_nodes=migrated_nodes+1
				comp_grid(xk,yk,zk)=comp(c8,4)
			else
				repeated_nodes=repeated_nodes+1
			end if
		end do
		
		write(3,*) "Number of samples located out of the migration grid:"
		write(3,*) comp_out_of_bounds
		write(3,*) "Number of migrated samples:"
		write(3,*) migrated_nodes
		write(3,*) "Number of repeated samples (not migrated and lost):"
		write(3,*) repeated_nodes
		
		do c9=1,ti_num
			do c10=1,blocks_z_mti
				do c11=1,blocks_y_mti
					do c12=1,blocks_x_mti
						mti(c12,c11,c10,c9)=-999
					end do
				end do
			end do
		end do
		
		
		do c13=1,blocks_z_mcomp_grid
			do c14=1,blocks_y_mcomp_grid
				do c15=1,blocks_x_mcomp_grid
					mcomp_grid(c15,c14,c13)=-999
				end do
			end do
		end do
		
		do c16=1,ti_num
			do c17=1,blocks_z_ti
				do c18=1,blocks_y_ti
					do c19=1,blocks_x_ti
						mti(c19+search_rad_x,c18+search_rad_y,c17+search_rad_z,c16)=ti(c19,c18,c17,c16)
					end do
				end do
			end do
		end do
		
		do c20=1,blocks_z_comp
			do c21=1,blocks_y_comp
				do c22=1,blocks_x_comp
					mcomp_grid(c22+search_rad_x,c21+search_rad_y,c20+search_rad_z)=comp_grid(c22,c21,c20)
				end do
			end do
		end do
		
		aux1=1
		
		do c23=1,zt
			do c24=1,yt
				do c25=1,xt
					tmp(aux1,1)=c25-search_rad_x-1
					tmp(aux1,2)=c24-search_rad_y-1
					tmp(aux1,3)=c23-search_rad_z-1
					tmp(aux1,4)=0
					tmp(aux1,5)=sqrt((tmp(aux1,1)*tmp(aux1,1))+(tmp(aux1,2)*tmp(aux1,2))+(tmp(aux1,3)*tmp(aux1,3)))
					aux1=aux1+1
				end do
			end do
		end do
		
		aux2=0
		orderer=1000000
		
		do c41=1,nt
			do c42=1,nt
				if (tmp(c42,5)<orderer) then
					orderer=tmp(c42,5)
					aux2=c42
					tmp_ord(c41,1)=tmp(c42,1)
					tmp_ord(c41,2)=tmp(c42,2)
					tmp_ord(c41,3)=tmp(c42,3)
					tmp_ord(c41,4)=tmp(c42,4)
				end if
			end do
			tmp(aux2,5)=1100000
			orderer=1000000
			aux2=0
		end do
				
				
				
		aux1=1
		
		do c26=1,max_conditioning_nodes
			c_event(c26,1)=-111
			c_event(c26,2)=-111
			c_event(c26,3)=-111
			c_event(c26,4)=-111
		end do
		
		startx=1+search_rad_x
		starty=1+search_rad_y
		startz=1+search_rad_z
		visited_node=0
		cond_node=0
		valid_events=0
		invalid_events=0
		total_events=0
		checker=0
		
		do c27=1,blocks_z_mcomp_grid
			do c28=1,blocks_y_mcomp_grid
				do c29=1,blocks_x_mcomp_grid
					if (mcomp_grid(c29,c28,c27)==-999) then
						cycle
					else
						do c30=1,nt
							offx=tmp_ord(c30,1)
							offy=tmp_ord(c30,2)
							offz=tmp_ord(c30,3)
							coordtempx=c29+offx
							coordtempy=c28+offy
							coordtempz=c27+offz
							if (tmp_ord(c30,4)==1) then
								cycle
							else if (tmp_ord(c30,4)==0 .and. mcomp_grid(coordtempx,coordtempy,coordtempz)==-111) then
								visited_node=visited_node+1
								tmp_ord(c30,4)=1
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)==-999) then
								cycle
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)/=-999) then
								visited_node=visited_node+1
								cond_node=cond_node+1
								tmp_ord(c30,4)=1
								c_event(cond_node,1)=offx
								c_event(cond_node,2)=offy
								c_event(cond_node,3)=offz
								c_event(cond_node,4)=mcomp_grid(coordtempx,coordtempy,coordtempz)
							end if
							if (cond_node==max_conditioning_nodes) then
								exit
							end if
						end do
						checker=0
						auxr2=cond_event_tolerance*cond_node
						min_equal_nodes=cond_node-floor(auxr2)
						if (cond_node<min_conditioning_nodes) then
							invalid_events=invalid_events+1
						else
							valid_events=valid_events+1
							do c33=1,ti_num
								do c34=startz,blocks_z_mti
									do c35=starty,blocks_y_mti
										do c36=startx,blocks_x_mti
											do c37=1,cond_node
												tioffx=c_event(c37,1)
												tioffy=c_event(c37,2)
												tioffz=c_event(c37,3)
												tiofffacie=c_event(c37,4)
												coordtix=c36+tioffx
												coordtiy=c35+tioffy
												coordtiz=c34+tioffz
												if (var_type==0) then
													if (mti(coordtix,coordtiy,coordtiz,c33)==tiofffacie) then
														equal_nodes=equal_nodes+1
													end if
												else if (var_type==1) then
													if (mti(coordtix,coordtiy,coordtiz,c33)>=(tiofffacie-acc_thresh).and.&
													mti(coordtix,coordtiy,coordtiz,c33)<=(tiofffacie+acc_thresh)) then
														equal_nodes=equal_nodes+1
													end if
												end if
											end do
											if (equal_nodes>=min_equal_nodes) then
												cm(1,c33)=cm(1,c33)+1
												checker=checker+1
											end if
											equal_nodes=0
										end do
									end do
								end do		
							end do
							if (checker>0) then
								do c45=1,ti_num
									aux4=aux4+cm(1,c45)
								end do
								do c46=1,ti_num
									cm2(1,c46)=cm2(1,c46)+(cm(1,c46)/aux4)
								end do
								aux4=0
								do c47=1,ti_num
									cm(1,c47)=0
								end do
							end if
						end if
					end if
					visited_node=0
					cond_node=0					
					do c31=1,nt
						if (c31==mid_n_tmp) then
							tmp_ord(c31,4)=0
						else
							tmp_ord(c31,4)=0
						end if
					end do
					
					do c32=1,max_conditioning_nodes
						c_event(c32,1)=-111
						c_event(c32,2)=-111
						c_event(c32,3)=-111
						c_event(c32,4)=-111
					end do
					
					total_events=valid_events+invalid_events
					write(*,*) total_events
				
				end do
			end do
		end do
		
		write(3,*) "Valid conditioning events:"
		write(3,*) valid_events
		write(3,*) "Invalid conditioning events:"
		write(3,*) invalid_events
		
		aux3=0
		
		do c38=1,ti_num
			aux3=aux3+cm2(1,c38)
		end do
		
		do c39=1,ti_num
			ncm(1,c39)=cm2(1,c39)/aux3
		end do
		
		write(3,*) "Relative compatibility with training images (exhaustive scanning):"
		write(3,*) (ncm(1,c51),c51=1,ti_num)			
		
		close(1)
		close(2)
		close(3)
	end if
	
		
	if (mode==0.and.ti_check_proportion>0) then
		comp_out_of_bounds=0
		migrated_nodes=0
		repeated_nodes=0
		
		
		blocks_in_ti=blocks_x_ti*blocks_y_ti*blocks_z_ti
		blocks_in_comp=blocks_x_comp*blocks_y_comp*blocks_z_comp
		image_min_x=image_min_x_centroid-(ti_x_block_size/2)
		image_min_y=image_min_y_centroid-(ti_y_block_size/2)
		image_min_z=image_min_z_centroid-(ti_z_block_size/2)
		comp_min_x=comp_min_x_centroid-(comp_x_block_size/2)
		comp_min_y=comp_min_y_centroid-(comp_y_block_size/2)
		comp_min_z=comp_min_z_centroid-(comp_z_block_size/2)
		image_max_x=image_min_x_centroid+(ti_x_block_size*(blocks_x_ti-1))+(ti_x_block_size/2)
		image_max_y=image_min_y_centroid+(ti_y_block_size*(blocks_y_ti-1))+(ti_y_block_size/2)
		image_max_z=image_min_z_centroid+(ti_z_block_size*(blocks_z_ti-1))+(ti_z_block_size/2)
		comp_max_x=comp_min_x_centroid+(comp_x_block_size*(blocks_x_comp-1))+(comp_x_block_size/2)
		comp_max_y=comp_min_y_centroid+(comp_y_block_size*(blocks_y_comp-1))+(comp_y_block_size/2)
		comp_max_z=comp_min_z_centroid+(comp_z_block_size*(blocks_z_comp-1))+(comp_z_block_size/2)
		blocks_x_mti=blocks_x_ti+2*search_rad_x
		blocks_y_mti=blocks_y_ti+2*search_rad_y
		blocks_z_mti=blocks_z_ti+2*search_rad_z
		blocks_in_mti=blocks_x_mti*blocks_y_mti*blocks_z_mti
		blocks_x_mcomp_grid=blocks_x_comp+2*search_rad_x
		blocks_y_mcomp_grid=blocks_y_comp+2*search_rad_y
		blocks_z_mcomp_grid=blocks_z_comp+2*search_rad_z
		blocks_in_mcomp_grid=blocks_x_mcomp_grid*blocks_y_mcomp_grid*blocks_z_mcomp_grid
		xt=search_rad_x*2+1
		yt=search_rad_y*2+1
		zt=search_rad_z*2+1
		nt=xt*yt*zt
		auxr1=nt/2
		mid_n_tmp=floor(auxr1)+1
		blocks_to_check=blocks_in_ti*ti_check_proportion
		
		
		call init_random_seed()
		
		do seed_counter=1,random_seed
			call random_number(random)
		end do
		
		open(1,file=comp_file)
		read(1,*,err=99)
		read(1,*,err=99) nvar
		do read_file_counter=1,nvar
			read(1,*,err=99)
		end do
		
		open(2,file=image_file)
		read(2,*,err=97)
		read(2,*,err=97) ti_num
		do read_file_counter=1,ti_num
			read(2,*,err=97)
		end do
		
		allocate (ti(blocks_x_ti,blocks_y_ti,blocks_z_ti,ti_num))
		allocate (tif(ti_num))
		allocate (mti(blocks_x_mti,blocks_y_mti,blocks_z_mti,ti_num))
		allocate (comp_grid(blocks_x_comp,blocks_y_comp,blocks_z_comp))
		allocate (mcomp_grid(blocks_x_mcomp_grid,blocks_y_mcomp_grid,blocks_z_mcomp_grid))
		allocate (tmp(nt,5))
		allocate (tmp_ord(nt,4))
		allocate (cm(1,ti_num))
		allocate (ncm(1,ti_num))
		allocate (cm2(1,ti_num))
		
		do c40=1,ti_num
			cm(1,c40)=0
			ncm(1,c40)=0
		end do


		open(3,file=output_file)
		
		do c1=1,comp_num_in_file
			read(1,*,err=99) xcomp,ycomp,zcomp,fcomp
			comp(c1,1)=xcomp
			comp(c1,2)=ycomp
			comp(c1,3)=zcomp
			comp(c1,4)=fcomp
		end do
		
		do c2=1,blocks_z_ti
			do c3=1,blocks_y_ti
				do c4=1,blocks_x_ti
					read(2,*,err=97) (tif(c49),c49=1,ti_num)
					do c50=1,ti_num
						ti(c4,c3,c2,c50)=tif(c50)
					end do
				end do
			end do
		end do
		
		do c5=1,blocks_z_comp
			do c6=1,blocks_y_comp
				do c7=1,blocks_x_comp
					comp_grid(c7,c6,c5)=-111
				end do
			end do
		end do
		
		do c8=1,comp_num_in_file
			if (comp(c8,1)<comp_min_x .or. comp(c8,1)>comp_max_x .or. comp(c8,2)<comp_min_y .or. comp(c8,2)>comp_max_y) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			if (comp(c8,3)<comp_min_z .or. comp(c8,3)>comp_max_z) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			xi=comp(c8,1)-comp_min_x
			yi=comp(c8,2)-comp_min_y
			zi=comp(c8,3)-comp_min_z
			xj=xi/comp_x_block_size+1
			yj=yi/comp_y_block_size+1
			zj=zi/comp_z_block_size+1
			xk=floor(xj)
			yk=floor(yj)
			zk=floor(zj)
			if (comp_grid(xk,yk,zk)==-111) then
				migrated_nodes=migrated_nodes+1
				comp_grid(xk,yk,zk)=comp(c8,4)
			else
				repeated_nodes=repeated_nodes+1
			end if
		end do
		
		write(3,*) "Number of samples located out of the training image grid:"
		write(3,*) comp_out_of_bounds
		write(3,*) "Number of migrated samples:"
		write(3,*) migrated_nodes
		write(3,*) "Number of repeated samples (not migrated and lost):"
		write(3,*) repeated_nodes
		
		do c9=1,ti_num
			do c10=1,blocks_z_mti
				do c11=1,blocks_y_mti
					do c12=1,blocks_x_mti
						mti(c12,c11,c10,c9)=-999
					end do
				end do
			end do
		end do
		
		
		do c13=1,blocks_z_mcomp_grid
			do c14=1,blocks_y_mcomp_grid
				do c15=1,blocks_x_mcomp_grid
					mcomp_grid(c15,c14,c13)=-999
				end do
			end do
		end do
		
		do c16=1,ti_num
			do c17=1,blocks_z_ti
				do c18=1,blocks_y_ti
					do c19=1,blocks_x_ti
						mti(c19+search_rad_x,c18+search_rad_y,c17+search_rad_z,c16)=ti(c19,c18,c17,c16)
					end do
				end do
			end do
		end do
		
		do c20=1,blocks_z_ti
			do c21=1,blocks_y_ti
				do c22=1,blocks_x_ti
					mcomp_grid(c22+search_rad_x,c21+search_rad_y,c20+search_rad_z)=comp_grid(c22,c21,c20)
				end do
			end do
		end do
		
		aux1=1
		
		do c23=1,zt
			do c24=1,yt
				do c25=1,xt
					tmp(aux1,1)=c25-search_rad_x-1
					tmp(aux1,2)=c24-search_rad_y-1
					tmp(aux1,3)=c23-search_rad_z-1
					tmp(aux1,4)=0
					tmp(aux1,5)=sqrt((tmp(aux1,1)*tmp(aux1,1))+(tmp(aux1,2)*tmp(aux1,2))+(tmp(aux1,3)*tmp(aux1,3)))
					aux1=aux1+1
				end do
			end do
		end do
		
		aux2=0
		orderer=1000000
		
		do c41=1,nt
			do c42=1,nt
				if (tmp(c42,5)<orderer) then
					orderer=tmp(c42,5)
					aux2=c42
					tmp_ord(c41,1)=tmp(c42,1)
					tmp_ord(c41,2)=tmp(c42,2)
					tmp_ord(c41,3)=tmp(c42,3)
					tmp_ord(c41,4)=tmp(c42,4)
				end if
			end do
			tmp(aux2,5)=1100000
			orderer=1000000
			aux2=0
		end do			
				
				
		aux1=1
		
		do c26=1,max_conditioning_nodes
			c_event(c26,1)=-111
			c_event(c26,2)=-111
			c_event(c26,3)=-111
			c_event(c26,4)=-111
		end do
		
		startx=1+search_rad_x
		starty=1+search_rad_y
		startz=1+search_rad_z
		visited_node=0
		cond_node=0
		valid_events=0
		invalid_events=0
		total_events=0
		checker=0
		
		do c27=1,blocks_z_mcomp_grid
			do c28=1,blocks_y_mcomp_grid
				do c29=1,blocks_x_mcomp_grid
					if (mcomp_grid(c29,c28,c27)==-999) then
						cycle
					else
						do c30=1,nt
							offx=tmp_ord(c30,1)
							offy=tmp_ord(c30,2)
							offz=tmp_ord(c30,3)
							coordtempx=c29+offx
							coordtempy=c28+offy
							coordtempz=c27+offz
							if (tmp_ord(c30,4)==1) then
								cycle
							else if (tmp_ord(c30,4)==0 .and. mcomp_grid(coordtempx,coordtempy,coordtempz)==-111) then
								visited_node=visited_node+1
								tmp_ord(c30,4)=1
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)==-999) then
								cycle
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)/=-999) then
								visited_node=visited_node+1
								cond_node=cond_node+1
								tmp_ord(c30,4)=1
								c_event(cond_node,1)=offx
								c_event(cond_node,2)=offy
								c_event(cond_node,3)=offz
								c_event(cond_node,4)=mcomp_grid(coordtempx,coordtempy,coordtempz)
							end if
							if (cond_node==max_conditioning_nodes) then
								exit
							end if
						end do
						auxr2=cond_event_tolerance*cond_node
						min_equal_nodes=cond_node-floor(auxr2)
						if (cond_node<min_conditioning_nodes) then
							invalid_events=invalid_events+1
						else
							valid_events=valid_events+1
							do c44=1,blocks_to_check
								call random_number(random)
								randomx=random*blocks_x_ti
								rintegerx=floor(randomx)+startx
								call random_number(random)
								randomy=random*blocks_y_ti
								rintegery=floor(randomy)+starty
								call random_number(random)
								randomz=random*blocks_z_ti
								rintegerz=floor(randomz)+startz
								do c33=1,ti_num
									do c37=1,cond_node
										tioffx=c_event(c37,1)
										tioffy=c_event(c37,2)
										tioffz=c_event(c37,3)
										tiofffacie=c_event(c37,4)
										coordtix=rintegerx+tioffx
										coordtiy=rintegery+tioffy
										coordtiz=rintegerz+tioffz
										if (var_type==0) then
											if (mti(coordtix,coordtiy,coordtiz,c33)==tiofffacie) then
												equal_nodes=equal_nodes+1
											end if
										else if (var_type==1) then
											if (mti(coordtix,coordtiy,coordtiz,c33)>=(tiofffacie-acc_thresh).and.&
											mti(coordtix,coordtiy,coordtiz,c33)<=(tiofffacie+acc_thresh)) then
												equal_nodes=equal_nodes+1
											end if
										end if
									end do
									if (equal_nodes>=min_equal_nodes) then
										cm(1,c33)=cm(1,c33)+1
										checker=1
									end if
									equal_nodes=0
								end do
								if (checker==1) then
									checker=0
									exit
								end if
							end do
						
						end if
					end if
					visited_node=0
					cond_node=0					
					do c31=1,nt
						if (c31==mid_n_tmp) then
							tmp_ord(c31,4)=0
						else
							tmp_ord(c31,4)=0
						end if
					end do
					do c32=1,max_conditioning_nodes
						c_event(c32,1)=-111
						c_event(c32,2)=-111
						c_event(c32,3)=-111
						c_event(c32,4)=-111
					end do
					
					total_events=valid_events+invalid_events
					write(*,*) total_events
				
				end do
			end do
		end do
		
		write(3,*) "Valid conditioning events:"
		write(3,*) valid_events
		write(3,*) "Invalid conditioning events:"
		write(3,*) invalid_events
		
		aux3=0
		
		do c38=1,ti_num
			aux3=aux3+cm(1,c38)
		end do
		
		do c39=1,ti_num
			ncm(1,c39)=cm(1,c39)/aux3
		end do
		
		write(3,*) "Relative compatibility with training images (direct sampling):"
		write(3,*) (ncm(1,c51),c51=1,ti_num)	
		write(3,*) "Occurrences of conditioning events in each training image:"
		write(3,*) (cm(1,c52),c52=1,ti_num)
			
		close(1)
		close(2)
		close(3)
	
	end if
	
	if (mode==1) then
		comp_out_of_bounds=0
		migrated_nodes=0
		repeated_nodes=0
		
		
		blocks_in_ti=blocks_x_ti*blocks_y_ti*blocks_z_ti
		blocks_in_comp=blocks_x_comp*blocks_y_comp*blocks_z_comp
		image_min_x=image_min_x_centroid-(ti_x_block_size/2)
		image_min_y=image_min_y_centroid-(ti_y_block_size/2)
		image_min_z=image_min_z_centroid-(ti_z_block_size/2)
		comp_min_x=comp_min_x_centroid-(comp_x_block_size/2)
		comp_min_y=comp_min_y_centroid-(comp_y_block_size/2)
		comp_min_z=comp_min_z_centroid-(comp_z_block_size/2)
		image_max_x=image_min_x_centroid+(ti_x_block_size*(blocks_x_ti-1))+(ti_x_block_size/2)
		image_max_y=image_min_y_centroid+(ti_y_block_size*(blocks_y_ti-1))+(ti_y_block_size/2)
		image_max_z=image_min_z_centroid+(ti_z_block_size*(blocks_z_ti-1))+(ti_z_block_size/2)
		comp_max_x=comp_min_x_centroid+(comp_x_block_size*(blocks_x_comp-1))+(comp_x_block_size/2)
		comp_max_y=comp_min_y_centroid+(comp_y_block_size*(blocks_y_comp-1))+(comp_y_block_size/2)
		comp_max_z=comp_min_z_centroid+(comp_z_block_size*(blocks_z_comp-1))+(comp_z_block_size/2)
		blocks_x_mti=blocks_x_ti+2*search_rad_x
		blocks_y_mti=blocks_y_ti+2*search_rad_y
		blocks_z_mti=blocks_z_ti+2*search_rad_z
		blocks_in_mti=blocks_x_mti*blocks_y_mti*blocks_z_mti
		blocks_x_mcomp_grid=blocks_x_comp+2*search_rad_x
		blocks_y_mcomp_grid=blocks_y_comp+2*search_rad_y
		blocks_z_mcomp_grid=blocks_z_comp+2*search_rad_z
		blocks_in_mcomp_grid=blocks_x_mcomp_grid*blocks_y_mcomp_grid*blocks_z_mcomp_grid
		xt=search_rad_x*2+1
		yt=search_rad_y*2+1
		zt=search_rad_z*2+1
		nt=xt*yt*zt
		auxr1=nt/2
		mid_n_tmp=floor(auxr1)+1
		blocks_to_check=blocks_in_ti*ti_check_proportion
		
		
		call init_random_seed()
		
		do seed_counter=1,random_seed
			call random_number(random)
		end do
		
		open(1,file=comp_file)
		read(1,*,err=99)
		read(1,*,err=99) nvar
		do read_file_counter=1,nvar
			read(1,*,err=99)
		end do
		
		open(2,file=image_file)
		read(2,*,err=97)
		read(2,*,err=97) ti_num
		if (ti_num<column) then
			goto 96
		end if
		do read_file_counter=1,ti_num
			read(2,*,err=97)
		end do
		
		allocate (ti(blocks_x_ti,blocks_y_ti,blocks_z_ti,ti_num))
		allocate (tif(ti_num))
		allocate (mti(blocks_x_mti,blocks_y_mti,blocks_z_mti,ti_num))
		allocate (comp_grid(blocks_x_comp,blocks_y_comp,blocks_z_comp))
		allocate (mcomp_grid(blocks_x_mcomp_grid,blocks_y_mcomp_grid,blocks_z_mcomp_grid))
		allocate (tmp(nt,5))
		allocate (tmp_ord(nt,4))
		allocate (cm(1,ti_num))
		allocate (ncm(1,ti_num))
		allocate (cm2(1,ti_num))
		
		do c40=1,ti_num
			cm(1,c40)=0
			ncm(1,c40)=0
		end do


		open(3,file=output_file)
		
		do c1=1,comp_num_in_file
			read(1,*,err=99) xcomp,ycomp,zcomp,fcomp
			comp(c1,1)=xcomp
			comp(c1,2)=ycomp
			comp(c1,3)=zcomp
			comp(c1,4)=fcomp
		end do
		
		do c2=1,blocks_z_ti
			do c3=1,blocks_y_ti
				do c4=1,blocks_x_ti
					read(2,*,err=97) (tif(c49),c49=1,ti_num)
					do c50=1,ti_num
						ti(c4,c3,c2,c50)=tif(c50)
					end do
				end do
			end do
		end do
		
		do c5=1,blocks_z_comp
			do c6=1,blocks_y_comp
				do c7=1,blocks_x_comp
				comp_grid(c7,c6,c5)=-111
				end do
			end do
		end do
		
		do c8=1,comp_num_in_file
			if (comp(c8,1)<comp_min_x .or. comp(c8,1)>comp_max_x .or. comp(c8,2)<comp_min_y .or. comp(c8,2)>comp_max_y) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			if (comp(c8,3)<comp_min_z .or. comp(c8,3)>comp_max_z) then
				comp_out_of_bounds=comp_out_of_bounds+1
				cycle
			end if
			xi=comp(c8,1)-comp_min_x
			yi=comp(c8,2)-comp_min_y
			zi=comp(c8,3)-comp_min_z
			xj=xi/comp_x_block_size+1
			yj=yi/comp_y_block_size+1
			zj=zi/comp_z_block_size+1
			xk=floor(xj)
			yk=floor(yj)
			zk=floor(zj)
			if (comp_grid(xk,yk,zk)==-111) then
				migrated_nodes=migrated_nodes+1
				comp_grid(xk,yk,zk)=comp(c8,4)
			else
				repeated_nodes=repeated_nodes+1
			end if
		end do
		
		write(3,*) "Number of samples located out of the training image grid:"
		write(3,*) comp_out_of_bounds
		write(3,*) "Number of migrated samples:"
		write(3,*) migrated_nodes
		write(3,*) "Number of repeated samples (not migrated and lost):"
		write(3,*) repeated_nodes
		
		do c9=1,ti_num
			do c10=1,blocks_z_mti
				do c11=1,blocks_y_mti
					do c12=1,blocks_x_mti
						mti(c12,c11,c10,c9)=-999
					end do
				end do
			end do
		end do
		
		
		do c13=1,blocks_z_mcomp_grid
			do c14=1,blocks_y_mcomp_grid
				do c15=1,blocks_x_mcomp_grid
					mcomp_grid(c15,c14,c13)=-999
				end do
			end do
		end do
		
		do c17=1,blocks_z_ti
			do c18=1,blocks_y_ti
				do c19=1,blocks_x_ti
					mti(c19+search_rad_x,c18+search_rad_y,c17+search_rad_z,1)=ti(c19,c18,c17,column)
				end do
			end do
		end do
		
		do c20=1,blocks_z_ti
			do c21=1,blocks_y_ti
				do c22=1,blocks_x_ti
					mcomp_grid(c22+search_rad_x,c21+search_rad_y,c20+search_rad_z)=comp_grid(c22,c21,c20)
				end do
			end do
		end do
		
		aux1=1
		
		do c23=1,zt
			do c24=1,yt
				do c25=1,xt
					tmp(aux1,1)=c25-search_rad_x-1
					tmp(aux1,2)=c24-search_rad_y-1
					tmp(aux1,3)=c23-search_rad_z-1
					tmp(aux1,4)=0
					tmp(aux1,5)=sqrt((tmp(aux1,1)*tmp(aux1,1))+(tmp(aux1,2)*tmp(aux1,2))+(tmp(aux1,3)*tmp(aux1,3)))
					aux1=aux1+1
				end do
			end do
		end do
		
		aux2=0
		orderer=1000000
		
		do c41=1,nt
			do c42=1,nt
				if (tmp(c42,5)<orderer) then
					orderer=tmp(c42,5)
					aux2=c42
					tmp_ord(c41,1)=tmp(c42,1)
					tmp_ord(c41,2)=tmp(c42,2)
					tmp_ord(c41,3)=tmp(c42,3)
					tmp_ord(c41,4)=tmp(c42,4)
				end if
			end do
			tmp(aux2,5)=1100000
			orderer=1000000
			aux2=0
		end do			
				
				
		aux1=1
		
		do c26=1,max_conditioning_nodes
			c_event(c26,1)=-111
			c_event(c26,2)=-111
			c_event(c26,3)=-111
			c_event(c26,4)=-111
		end do
		
		startx=1+search_rad_x
		starty=1+search_rad_y
		startz=1+search_rad_z
		visited_node=0
		cond_node=0
		valid_events=0
		invalid_events=0
		total_events=0
		checker=0
		
		do c27=1,blocks_z_mcomp_grid
			do c28=1,blocks_y_mcomp_grid
				do c29=1,blocks_x_mcomp_grid
					if (mcomp_grid(c29,c28,c27)==-999) then
						cycle
					else
						do c30=1,nt
							offx=tmp_ord(c30,1)
							offy=tmp_ord(c30,2)
							offz=tmp_ord(c30,3)
							coordtempx=c29+offx
							coordtempy=c28+offy
							coordtempz=c27+offz
							if (tmp_ord(c30,4)==1) then
								cycle
							else if (tmp_ord(c30,4)==0 .and. mcomp_grid(coordtempx,coordtempy,coordtempz)==-111) then
								visited_node=visited_node+1
								tmp_ord(c30,4)=1
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)==-999) then
								cycle
							else if (mcomp_grid(coordtempx,coordtempy,coordtempz)/=-999) then
								visited_node=visited_node+1
								cond_node=cond_node+1
								tmp_ord(c30,4)=1
								c_event(cond_node,1)=offx
								c_event(cond_node,2)=offy
								c_event(cond_node,3)=offz
								c_event(cond_node,4)=mcomp_grid(coordtempx,coordtempy,coordtempz)
							end if
							if (cond_node==max_conditioning_nodes) then
								exit
							end if
						end do
						auxr2=cond_event_tolerance*cond_node
						min_equal_nodes=cond_node-floor(auxr2)
						if (cond_node<min_conditioning_nodes) then
							invalid_events=invalid_events+1
						else
							valid_events=valid_events+1
							do c44=1,blocks_to_check
								call random_number(random)
								randomx=random*blocks_x_ti
								rintegerx=floor(randomx)+startx
								call random_number(random)
								randomy=random*blocks_y_ti
								rintegery=floor(randomy)+starty
								call random_number(random)
								randomz=random*blocks_z_ti
								rintegerz=floor(randomz)+startz
								do c37=1,cond_node
									tioffx=c_event(c37,1)
									tioffy=c_event(c37,2)
									tioffz=c_event(c37,3)
									tiofffacie=c_event(c37,4)
									coordtix=rintegerx+tioffx
									coordtiy=rintegery+tioffy
									coordtiz=rintegerz+tioffz
									if (var_type==0) then
										if (mti(coordtix,coordtiy,coordtiz,1)==tiofffacie) then
											equal_nodes=equal_nodes+1
										end if
									else if (var_type==1) then
										if (mti(coordtix,coordtiy,coordtiz,1)>=(tiofffacie-acc_thresh).and.&
										mti(coordtix,coordtiy,coordtiz,1)<=(tiofffacie+acc_thresh)) then
											equal_nodes=equal_nodes+1
										end if
									end if
								end do
								if (equal_nodes>=min_equal_nodes) then
									cm(1,1)=cm(1,1)+1
									checker=1
								end if
								equal_nodes=0
								if (checker==1) then
									checker=0
									exit
								end if
							end do
						
						end if
					end if
					visited_node=0
					cond_node=0					
					do c31=1,nt
						if (c31==mid_n_tmp) then
							tmp_ord(c31,4)=0
						else
							tmp_ord(c31,4)=0
						end if
					end do
					do c32=1,max_conditioning_nodes
						c_event(c32,1)=-111
						c_event(c32,2)=-111
						c_event(c32,3)=-111
						c_event(c32,4)=-111
					end do
					
					total_events=valid_events+invalid_events
					write(*,*) total_events
				
				end do
			end do
		end do
		
		write(3,*) "Valid conditioning events:"
		write(3,*) valid_events
		write(3,*) "Invalid conditioning events:"
		write(3,*) invalid_events
		
		write(3,*) "Absolute compatibility with training images:"
		write(3,*) cm(1,1)/valid_events
		
			
		write(3,*) "Matching conditioning events:"
		write(3,*) cm(1,1)
				
		close(1)
		close(2)
		close(3)
	end if
		
		
	stop
95      stop 'ERROR compatibility mode has to be set to 0 (relative) or 1 (absolute)'
96	stop 'ERROR column to analyze does not exist in images file'
97	stop 'ERROR in training images file somewhere'
98	stop 'ERROR in parameters somewhere'
99	stop 'ERROR in  composites datafile somewhere'
	stop
	end
	
	
	SUBROUTINE init_random_seed()
	            INTEGER :: i, n, clock
	            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
	          
	            CALL RANDOM_SEED(size = n)
	            ALLOCATE(seed(n))
	          
	            CALL SYSTEM_CLOCK(COUNT=clock)
	          
	            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	            CALL RANDOM_SEED(PUT = seed)
	          
	            DEALLOCATE(seed)
	END SUBROUTINE
