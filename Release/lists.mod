	  !  ]   k820309    Ö          15.0        OšÎV                                                                                                           
       utils.F90 LISTS                   @                               'H                    #P                                                                           	            &                                                             @                                'H                    #P                                                                            
            &                                                             @                               'H                    #P    .                                                                                   &                                                                     @                                'X                    #COUNT    #DELTA 	   #CAPACITY 
   #ITEMS                                                                                                               	                                                              
                                                                               H             #REAL_POINTER              &                                                             @                                'X                    #COUNT    #DELTA    #CAPACITY    #ITEMS                                                                                                                                                                                                                                                            H             #STRING_POINTER              &                                           #         @                                                       #L              D                                      X               #TLIST_REALARR    #         @                                                      #L              D                                      X               #TLIST_REALARR    #         @                                                      #TLIST_REALARR_ADD%SIZE    #L    #P                                                   SIZE           D @                                    X               #TLIST_REALARR              
 @                                                  	              &                                           #         @                                                      #L    #C              D                                      X               #TLIST_REALARR                                                           #         @                                                       #L    #I              D                                      X               #TLIST_REALARR              
                                             #         @                                                      #TLIST_REALARR_SAVEBINARY%SIZE     #L !   #FID "                                                   SIZE                                            !     X               #TLIST_REALARR              
                                  "           #         @                                   #                    #L $   #FID %             D @                               $     X               #TLIST_REALARR              
                                  %           #         @                                   &                    #L '   #I (             D                                 '     X               #TLIST_REALARR              
                                  (           #         @                                   )                   #TLIST_REALARR_CONFIDVAL%INT *   #TLIST_REALARR_CONFIDVAL%MAX +   #TLIST_REALARR_CONFIDVAL%PRESENT ,   #L -   #IX .   #LIMFRAC /   #IX1 0   #IX2 1   #LOWER 2   #UPPER 3                                             *     INT                                           +     MAX                                           ,     PRESENT                                            -     X               #TLIST_REALARR              
  @                               .                     
                                  /     	                
 @                               0                     
 @                               1                     F @                               2     	                 F @                               3     	       #         @                                 4                    #ARR 5   #LIN 6   #R 7   #INDEX 8          @  D @                               5             H             p          1     1                   #REAL_POINTER              
                                  6                     
                                  7                     
  @                               8           #         @                                  9                    #L :             D                                 :     X               #TSTRINGLIST    #         @                                  ;                    #L <             D @                               <     X               #TSTRINGLIST    #         @                                   =                   #TSTRINGLIST_SETFROMSTRING%TRIM >   #TSTRINGLIST_SETFROMSTRING%VERIFY ?   #TSTRINGLIST_SETFROMSTRING%LEN_TRIM @   #TSTRINGLIST_SETFROMSTRING%PRESENT A   #L B   #S C   #VALID_CHARS_IN D                                             >     TRIM                                           ?     VERIFY                                           @     LEN_TRIM                                           A     PRESENT           D @                               B     X               #TSTRINGLIST              
  @                              C                    1           
 @                              D                    1 #         @                                  E                   #TSTRINGLIST_ADD%LEN_TRIM F   #L G   #P H                                             F     LEN_TRIM           D @                               G     X               #TSTRINGLIST              
  @                              H                    1 #         @                                  I                    #L J   #C K             D                                 J     X               #TSTRINGLIST                                               K            $         @                                 L                          #TSTRINGLIST_ITEM%SIZE M   #L N   #I O                                                     M     SIZE                                            N     X               #TSTRINGLIST              
                                  O           #         @                                   P                    #L Q   #I R             D                                 Q     X               #TSTRINGLIST              
                                  R           %         @                                 S                          #TSTRINGLIST_INDEXOF%LEN_TRIM T   #TSTRINGLIST_INDEXOF%SIZE U   #L V   #S W                                             T     LEN_TRIM                                           U     SIZE                                            V     X               #TSTRINGLIST              
  @                              W                    1 #         @                                 X                    #ARR Y   #LIN Z   #R [   #INDEX \          @  D @                               Y             H             p          1     1                   #DOUBLE_POINTER              
                                  Z                     
                                  [                     
  @                               \                        fn#fn    ¸   W       REAL_POINTER         a   REAL_POINTER%P    Ŗ  W       DOUBLE_POINTER !   ú     a   DOUBLE_POINTER%P      W       STRING_POINTER !   å     a   STRING_POINTER%P             TLIST_REALARR $      H   a   TLIST_REALARR%COUNT $   H  H   a   TLIST_REALARR%DELTA '     H   a   TLIST_REALARR%CAPACITY $   Ø  Ļ   a   TLIST_REALARR%ITEMS    ~         TSTRINGLIST "   ũ  H   a   TSTRINGLIST%COUNT "   E  H   a   TSTRINGLIST%DELTA %     H   a   TSTRINGLIST%CAPACITY "   Õ  ¨   a   TSTRINGLIST%ITEMS #   }  O       TLIST_REALARR_INIT %   Ė  [   a   TLIST_REALARR_INIT%L $   '  O       TLIST_REALARR_CLEAR &   v  [   a   TLIST_REALARR_CLEAR%L "   Ņ  r       TLIST_REALARR_ADD '   C	  =      TLIST_REALARR_ADD%SIZE $   	  [   a   TLIST_REALARR_ADD%L $   Û	     a   TLIST_REALARR_ADD%P *   g
  V       TLIST_REALARR_SETCAPACITY ,   Ŋ
  [   a   TLIST_REALARR_SETCAPACITY%L ,     @   a   TLIST_REALARR_SETCAPACITY%C %   X  V       TLIST_REALARR_DELETE '   Ž  [   a   TLIST_REALARR_DELETE%L '   	  @   a   TLIST_REALARR_DELETE%I )   I  {       TLIST_REALARR_SAVEBINARY .   Ä  =      TLIST_REALARR_SAVEBINARY%SIZE +     [   a   TLIST_REALARR_SAVEBINARY%L -   \  @   a   TLIST_REALARR_SAVEBINARY%FID )     X       TLIST_REALARR_READBINARY +   ô  [   a   TLIST_REALARR_READBINARY%L -   O  @   a   TLIST_REALARR_READBINARY%FID #     V       TLIST_REALARR_THIN %   å  [   a   TLIST_REALARR_THIN%L %   @  @   a   TLIST_REALARR_THIN%I (     ķ       TLIST_REALARR_CONFIDVAL ,   s  <      TLIST_REALARR_CONFIDVAL%INT ,   ¯  <      TLIST_REALARR_CONFIDVAL%MAX 0   ë  @      TLIST_REALARR_CONFIDVAL%PRESENT *   +  [   a   TLIST_REALARR_CONFIDVAL%L +     @   a   TLIST_REALARR_CONFIDVAL%IX 0   Æ  @   a   TLIST_REALARR_CONFIDVAL%LIMFRAC ,     @   a   TLIST_REALARR_CONFIDVAL%IX1 ,   F  @   a   TLIST_REALARR_CONFIDVAL%IX2 .     @   a   TLIST_REALARR_CONFIDVAL%LOWER .   Æ  @   a   TLIST_REALARR_CONFIDVAL%UPPER "     l       QUICKSORTARR_REAL &   r     a   QUICKSORTARR_REAL%ARR &     @   a   QUICKSORTARR_REAL%LIN $   H  @   a   QUICKSORTARR_REAL%R (     @   a   QUICKSORTARR_REAL%INDEX !   Č  O       TSTRINGLIST_INIT #     Y   a   TSTRINGLIST_INIT%L "   p  O       TSTRINGLIST_CLEAR $   ŋ  Y   a   TSTRINGLIST_CLEAR%L *           TSTRINGLIST_SETFROMSTRING /     =      TSTRINGLIST_SETFROMSTRING%TRIM 1   X  ?      TSTRINGLIST_SETFROMSTRING%VERIFY 3     A      TSTRINGLIST_SETFROMSTRING%LEN_TRIM 2   Ø  @      TSTRINGLIST_SETFROMSTRING%PRESENT ,     Y   a   TSTRINGLIST_SETFROMSTRING%L ,   q  L   a   TSTRINGLIST_SETFROMSTRING%S 9   Ŋ  L   a   TSTRINGLIST_SETFROMSTRING%VALID_CHARS_IN     	  t       TSTRINGLIST_ADD )   }  A      TSTRINGLIST_ADD%LEN_TRIM "   ž  Y   a   TSTRINGLIST_ADD%L "     L   a   TSTRINGLIST_ADD%P (   c  V       TSTRINGLIST_SETCAPACITY *   š  Y   a   TSTRINGLIST_SETCAPACITY%L *     @   a   TSTRINGLIST_SETCAPACITY%C !   R         TSTRINGLIST_ITEM &   Ķ  =      TSTRINGLIST_ITEM%SIZE #     Y   a   TSTRINGLIST_ITEM%L #   i  @   a   TSTRINGLIST_ITEM%I #   Š  V       TSTRINGLIST_DELETE %   ˙  Y   a   TSTRINGLIST_DELETE%L %   X  @   a   TSTRINGLIST_DELETE%I $            TSTRINGLIST_INDEXOF -   6  A      TSTRINGLIST_INDEXOF%LEN_TRIM )   w  =      TSTRINGLIST_INDEXOF%SIZE &   ´  Y   a   TSTRINGLIST_INDEXOF%L &     L   a   TSTRINGLIST_INDEXOF%S    Y  l       QUICKSORTARR !   Å     a   QUICKSORTARR%ARR !   ]   @   a   QUICKSORTARR%LIN       @   a   QUICKSORTARR%R #   Ũ   @   a   QUICKSORTARR%INDEX 