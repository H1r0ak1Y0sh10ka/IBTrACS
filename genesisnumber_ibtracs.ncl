;**********************************************************
;IBTRACS: Storm Cyclogenisis 
;**********************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/Volumes/HDD01/research/ncl/wrf/contributed_beta_640.ncl"
;**********************************************************

;********************
; get data
;********************
  dir  = "/Volumes/RAID01/research/DATA/ibtracs/"
  ;dir   = "/Volumes/Yoshi-HDD01/research/DATA/analysis/ibtracs/"
  ;data = "Allstorms.ibtracs_all."
  data  = "Allstorms.ibtracs_wmo."
  ver  = "v03r09"
  f = addfile(dir+data+ver+".nc","r")

  START_YEAR = 1981
  END_YEAR   = 2015

  START_LAT = 0.0
;  LAT_DIS = 2.0
  END_LAT  = 45.0
  NUM_GROUP = 20


;  NUM_GROUP = floattoint((END_LAT-START_LAT)/LAT_DIS)
  tcgroup =new(NUM_GROUP,integer)
  tcgroup(:) = 0
  YYY =new(NUM_GROUP+1,float)
  name_group=new(NUM_GROUP,string)

;  LAT_MIN = 0
;  LAT_MAX = 40

  START_F    = coriolis_param(START_LAT)
  END_F      = coriolis_param(END_LAT)

  F_DIS = (END_F-START_F)/NUM_GROUP

  PI    = 3.14159265358979    ;3.14159          ; 3.14159265358979
  OMEGA = 7.2921e-5 

  do ggg=0,NUM_GROUP
    DAM=F_DIS*ggg+START_LAT
    PHI=asin(DAM/(2*OMEGA))
    YYY(ggg)=PHI*(180/PI)    ;round(PHI*(180/PI),3)

    if(ggg.gt.0.and.ggg.le.NUM_GROUP)then
     name_group(ggg-1)=round(YYY(ggg-1),3)+"-"+round(YYY(ggg),3)
    end if
   end do

   print(name_group)

;**********************************************************

  lat    = short2flt(f->lat_wmo)       ; short2flt(f->lat_for_mapping)
  lon    = short2flt(f->lon_wmo)       ;short2flt(f->lon_for_mapping)
  wind   = short2flt(f->wind_wmo)
  status = f->nature_wmo
  basin  = f->genesis_basin
  time   = f->season
  
  MAX_TIME = dimsizes(lat(0,:))
  NUM_STORM = dimsizes(lat(:,0))

  genesis_lat = new(NUM_STORM,float)
  genesis_lon = new(NUM_STORM,float)
  genesis_lat(:) = lat@_FillValue

  yyyy = (get_unique_values(time))
  nyyy = (dimsizes(yyyy))

  do nnn = 0,(NUM_STORM-1)
   chck=0
 
    do ttt = 0,(MAX_TIME-1)
;     print("STORM_Number:"+nnn+",TIME_STEP:"+ttt+",STATUS:"+status(nnn,ttt))
         xxx01=status(nnn,ttt)

       if(.not.ismissing(xxx01))then
          if(xxx01.eq.0.and.chck.eq.0)
          ; print("Genesis")
           ;  print("STORM_Number:"+nnn+",TIME_STEP:"+ttt+",STATUS:"+xxx01)
             genesis_lat(nnn) = lat(nnn,ttt) 
           chck=chck+1
          end if

       else
      ;   print("missing") 
       end if 

     end do

   end do


  do nnn = 0,NUM_STORM-1

   if(time(nnn).ge.START_YEAR.and.time(nnn).le.END_YEAR)then

   ;print(nnn+", This TCG lat is "+lat(nnn,0)+".")

   do ngg = 0,NUM_GROUP-1
    
     AAA = YYY(ngg)
     BBB = YYY(ngg+1)
  ;   print("From "+AAA+" to "+BBB+"." )
     xxx02 = genesis_lat(nnn)

   if(.not.ismissing(xxx02))then
    if(xxx02.ge.AAA.and.xxx02.lt.BBB)then
;    print("Genesis of "+nnn+": from "+AAA+" to "+BBB+"." )
    tcgroup(ngg)=tcgroup(ngg)+1
    end if
    end if

   end do ; ngg loop

   end if

  end do ; nnn loop

  print(tcgroup)
  MAX_TC = round((int2flt(max(tcgroup))*0.01),3)*100

;*********************************
; create plot
;*********************************
;---Using "newps" gives us better looking results with the grid lines.
  wks  = gsn_open_wks ("eps", "genesisnumber_int-"+NUM_GROUP+"groups")           ; open workstation

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

        res02 = True
        res02@gsnDraw      = False                    ; Turn off draw and frame so
        res02@gsnFrame     = False                    ; we can attach some text.
        res02@gsnScale = True      
        res02@gsnXYBarChart = True            ; Create bar plot

        res02@gsnXYBarChartColors = (/"red"/)	
        res02@tiYAxisString = "Number of TCG" ; y-axis label      
        res02@tiXAxisString = "Latitude" ; x-axis label           
        res02@tiMainString = START_YEAR+"-"+END_YEAR+"_"+F_DIS


        res02@trYMinF   =   0        ; min value on y-axis
        res02@trYMaxF   =  MAX_TC+100        ; max value on y-axis
        res02@trXMinF   =   -1        ; min value on y-axis
        res02@trXMaxF   =  (NUM_GROUP)        ; max value on y-axis

        res02@tmYLOn              = True
        res02@tmYLMode        = "Manual"
        res02@tmYMajorGrid                = True    ; Turn on grid lines
        res02@tmYLMinorPerMajor = 1
	    res02@tmYLTickStartF  = 0
        res02@tmYLTickEndF    = MAX_TC+100
        res02@tmYLTickSpacingF= 100
        res02@tmXBMode   = "Explicit"
        x = ispan(0,NUM_GROUP-1,1)
        res02@tmXBValues = x             ; choose first 13 timesteps
        res02@tmXBLabels = name_group
        res02@tmXBLabelAngleF = 45	; tilt the XB labels 45 degrees


; this controls the size and location of the second plot
       res02@gsnMaximize = True
       res02@vpWidthF      = 0.8
       res02@vpHeightF     = 0.2

;---------- First plot --------------------
        bar = gsn_csm_xy (wks,x,tcgroup,res02)

;
; Remove this text string, and create some new ones to add.
;
  txres                       = True
  amres                       = True
  txres@txBackgroundFillColor = "Black"
  txres@txFont                = "helvetica-bold"
  txres@txFontColor           = "white"
  txres@txFontHeightF         = 0.02
  txid1 = gsn_create_text(wks,"ALL:"+sum(tcgroup), txres)

  amres@amParallelPosF   = 0.5    ; This is the right edge of the plot.
  amres@amOrthogonalPosF = -0.5    ; This is the bottom edge of the plot.
;
; By default, the center of the string is what's placed at the position
; indicated by amParallelPosF and amOrthogonalPosF. You can use amJust
; to change this to any one of 9 positions: "CenterCenter" (default),
; "TopCenter", "TopRight", "CenterRight", "BottomRight", "BottomCenter",
; "BottomLeft", "CenterLeft", "TopLeft".
;
  amres@amJust           = "TopRight"
  annoid1 = gsn_add_annotation(bar, txid1, amres)

;**********************************************************
; add text labels
;**********************************************************  
  delete(txres)
  txres                     = True                 ; text mods desired
  txres@txFontHeightF       = 0.018                ; default size is HUGE!
  txres@txAngleF            = 52.                  ; text angle
  txres@txJust              = "CenterLeft"         ; puts text on top of bars
  
  do n = 0, NUM_GROUP-1
    gsn_text(wks,bar,tostring(tcgroup(n)),n,tcgroup(n)+50,txres) ; add labels
  end do

  draw(bar)    ; Both maps will get drawn.
  frame(wks)

cputime_min=get_cpu_time/60
print(sprintf("cputime=%f min",cputime_min))

exit
end


