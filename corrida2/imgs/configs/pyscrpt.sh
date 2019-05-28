for i in 1_5 2_0 2_5 3_0 3_5 4_0
do 
	for j in 1_0 1_5 2_0 2_5 3_0 3_5 4_0
	do 
		for k in 0_1 0_2 0_3 0_4 0_5 0_6 0_7 0_8 0_9 1_0
		do 
			python config.py l${i}r${j}rho${k}.xyz $1 
		done
	done
done
