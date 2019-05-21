path=`pwd`
echo $path
 for i in 1 1 
 do 
	 for j in 1 2 3
	 do
		 python time.py $i$j &
	 done
	 wait
 done
