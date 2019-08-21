all: gllvm2.so
	# Comment here preserves the prior tab
gllvm2.so: gllvm2.cpp
	Rscript --vanilla -e "TMB::compile('gllvm2.cpp','-O0 -g')"
  
clean:
	rm -rf *o
