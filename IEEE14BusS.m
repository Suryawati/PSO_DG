function [busdata,linedata,GenRestric] = IEEE14BusS

%        Bus Bus  Voltage Angle      ---Load----     ----------Generator---------   Static Mvar
%        No  code Mag.    Degree     MW    Mvar      MW         Mvar     Qmin     Qmax    +Qc/-Ql
busdata=[1	1   1.06     0          0       0	    232.4      -16.9	   0	   0	     0
        2	2	1.045	-4.98	    21.7    12.7	40	        42.4	  -40	  50	     0
        3	2	1.01	-12.72	    94.2    19	     0	        23.4	   0	  40	     0
        4	0	1.019	-10.33	    47.8    -3.9	 0	         0	       0	   0    	 0
        5	0	1.02	-8.78	     7.6     1.6	 0	         0	       0	   0	     0
        6	2	1.07	-14.22	    11.2	 7.5	 0	        12.2	  -6       24        0
        7	0	1.062	-13.37       0       0	     0	         0	       0	   0	     0
        8	2	1.09	-13.36       0       0	     0	        17.4	  -6       24        0
        9	0	1.056	-14.94	    29.5	16.6	 0	         0	       0	   0	     0.19
        10	0	1.051	-15.1        9	     5.8	 0	         0	       0	   0	     0
        11	0	1.057	-14.79	     3.5	 1.8	 0	         0	       0	   0	     0
        12	0	1.055	-15.07	     6.1	 1.6	 0	         0	       0	   0	     0
        13	0	1.05	-15.16	    13.5	 5.8	 0	         0	       0	   0	     0
        14	0	1.036	-16.04	    14.9	 5	     0	         0	       0	   0	     0];
        
%                                           Line code
%         Bus bus   R        X         1/2 B     = 1 for lines
%         nl  nr  p.u.      p.u.       p.u.     > 1 or < 1 tr. tap at bus nl
linedata=[1   2   0.01938   0.05917     0.0528          1
          1   5   0.05403   0.22304     0.0492          1
          2   3   0.04699   0.19797     0.0438          1
          2   4   0.05811   0.17632     0.0340          1
          2   5   0.05695   0.17388     0.0346          1
          3   4   0.06701   0.17103     0.0128          1
          4   5   0.01335   0.04211     0.0             1
          4   7   0.0       0.20912     0.0             1
          4   9   0.0       0.55618     0.0             1
          5   6   0.0       0.25202     0.0             1
          6  11   0.09498   0.19890     0.0             1
          6  12   0.12291   0.25581     0.0             1
          6  13   0.06615   0.13027     0.0             1
          7   8   0.0       0.17615     0.0             1
          7   9   0.0       0.11001     0.0             1
          9  10   0.03181   0.08450     0.0             1
          9  14   0.12711   0.27038     0.0             1
         10  11   0.08205   0.19207     0.0             1
         12  13   0.22092   0.19988     0.0             1
         13  14   0.17093   0.34802     0.0             1];
     
     GenRestric=[1 2 3 4 5 6 7 8];  %Bus tidak boleh di tempati DG untuk 14 Bus
     %TLoss=loadflow(busdata,linedata,1);