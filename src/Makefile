all: gllvm2.so gllvm3.so
	# Comment here preserves the prior tab
gllvm2.so: gllvm2.cpp
	Rscript --vanilla -e "TMB::compile('gllvm2.cpp','-O0 -g')"
	# Comment here preserves the prior tab
gllvm3.so: gllvm3.cpp
	Rscript --vanilla -e "TMB::compile('gllvm3.cpp','-O0 -g')"
	
clean:
	rm -rf *o
