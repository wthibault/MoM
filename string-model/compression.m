function gain = compression ( rmsLevel )
  compressionThreshold = -10
  compressionRatio = 10
	 if ( rmsLevel >  compressionThreshold )
	    gain = 1.0 - decibelsToLinear((rmsLevel-compressionThreshold)/compressionRatio)
         else
	    gain = 1.0
         endif
endfunction

  
