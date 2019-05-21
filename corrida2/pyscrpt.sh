path=$(dirname `pwd`) # la dirección del fichero arriba
for l in 1_5 2_0 2_5 3_0 3_5 4_0
do 
	cd lamb$l
	for i in 1_0 1_5 2_0 2_5 3_0 3_5 4_0
	do 
		cd r$i
		for j in 0_1 0_2 0_3 0_4 0_5 0_6 0_7 0_8 0_9 1_0
		do 
			cd rho00_$j
		        #cp $path/pyscrpts/config.py config.py
			#cp $path/pyscrpts/funcs.py funcs.py
			#cp $path/pyscrpts/qtree.py qtree.py
			#python config.py l${l}r${i}rho${j}.xyz $path/
			#rm funcs.py
			#rm qtree.py
			#rm config.py
			## rm lyr[0-9]*
			cp l${l}r${i}rho${j}.xyz ${path}/corrida2/imgs/.
			cd ..
		done
		cd ..
	done
	cd ..
done
