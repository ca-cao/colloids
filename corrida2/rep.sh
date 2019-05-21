path=`pwd` # la dirección donde están:
for l in 1_5 2_0 2_5 3_0 3_5 4_0
do 
	cd lamb$l
	for i in 1_0 1_5 2_0 2_5 3_0 3_5 4_0
	do 
	#if [ ! -d r$i ]
	#then	
	#	mkdir r$i 
	#fi
		cd r$i
		for j in 0_1 0_2 0_3 0_4 0_5 0_6 0_7 0_8 0_9 1_0
		do 
		#if [ ! -d rho00_$j ]
		#then 
		#	mkdir rho00_$j
		#fi
			cd rho00_$j
			#if [ ! -f data18.dat ]
			#then
				#rm  l${l}r${i}rho${j}.xyz
		        	#cp $path/t.f90 t.f90
				cp $path/input.in input.in
				#cp $path/random_seed.in random_seed.in
				#sed -i "19s!x!`echo $i | sed -e 's:_:.:'`!" t.f90  #cambia los parametros
				#sed -i "20s!x!`echo $j | sed -e 's:_:.:'`!" t.f90  #cambia los parametros
				#sed -i "27s!x!`echo $l | sed -e 's:_:.:'`!" t.f90  #cambia los parametros
				#sed -i "28s!x!`echo $l | sed -e 's:_:.:'`!" t.f90  #cambia los parametros
				#sed -i "37s!x!l${l}r${i}rho${j}.xyz!" t.f90  #cambia el nombre del archivo de salida
				#mv *.f90 l${l}r${i}rho${j}.f90
				#rm *.o
				#rm *.xyz
				gfortran l${l}r${i}rho${j}.f90 -o l${l}r${i}rho${j}.exe
				#./l${l}r${i}rho${j}.exe &
			#fi
			cd ..
		done
		wait
		cd ..
	done
	cd ..
done
