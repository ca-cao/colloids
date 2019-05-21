for i in $( ls -d lamb*/ ) 
do 
	cd  ${i}finalconfigs/
	for x in $( ls l* )
	do
		~/ovito-2.9.0-x86_64/bin/ovitos test_image.py $x
	done
	cd ../..
done
