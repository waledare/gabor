dyn.load("foo.so")

I <- complex(real=0, imaginary=1)
xs <- rnorm(5) + rnorm(5)*I
print(xs)

xls <- .C("foo",
          n=as.integer(5), 
          realxs=as.double(Re(xs)),
          imagxs=as.double(Im(xs)))
print(xs)
