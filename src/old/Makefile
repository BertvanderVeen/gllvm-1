all: qgllvm.so
	# Comment here preserves the prior tab
qgllvm.so: qgllvm.cpp
	Rscript --vanilla -e "TMB::compile('qgllvm.cpp','-O0 -g')"
	
clean:
	rm -rf *o
