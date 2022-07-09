testbuild:
	runhaskell Setup configure --user --enable-tests --enable-benchmarks
	runhaskell Setup build
	runhaskell Setup haddock
	runhaskell Setup test --show-details=streaming
