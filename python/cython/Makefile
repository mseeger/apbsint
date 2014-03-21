all:
	python setup.py build_ext --inplace 
	mv eptools_ext.so ../apbsint/.
	mv apbtest_ext.so ../apbsint/.
	mv ptannotate_ext.so ../apbsint/.

wa:
	python setup.py build_ext --inplace --workaround
	mv eptools_ext.so ../apbsint/.
	mv apbtest_ext.so ../apbsint/.
	mv apbtest_workaround_ext.so ../apbsint/.
	mv ptannotate_ext.so ../apbsint/.

clean:
	rm -rf build
	rm *.cpp
