function nlmap,x
  y=1.1441003e-11
  y=y*x-1.1310333e-08
  y=y*x+4.2562775e-06
  y=y*x-0.00077292738
  y=y*x+0.34734467
  y=y*x+0.18928926
  return,y
end

function map,f,xp,yp
  csz=size(xp)
  nl=csz[3]
  cml=total(f[0:99,*],1)/100.0
  dim=size(f)
  for x=0,dim[1]/2-1 do f[x,*]-=cml
  cmh=total(f[dim[1]-100:dim[1]-1,*],1)/100.0
  for x=dim[1]/2,dim[1]-1 do f[x,*]-=cmh
  f=nlmap(f)
  cube=dblarr(csz[1],csz[2],nl)
  for il=0,nl-1 do
cube[*,*,il]=interpolate(f,xp[*,*,il],yp[*,*,il],cubic=-0.5)
  return,reverse(transpose(reverse(cube,2),[1,0,2]),2)
end

function mkavg,imn,imx

 
avg=double(f0('/SST/2017.07.10/obs07/mihi_sp/ff3/image.'+string(imn,format='(I09)')+'.f0'))
  nn=1L
  for i=imn+3L,imx,3L do begin
  
avg+=double(f0('../../obs07/mihi_sp/ff3/image.'+string(i,format='(I09)')+'.f0'))
    nn+=1
  end
  return,avg/double(nn)
end

pro decmp
  xp=f0('xpos.f0')+1.0
  yp=f0('ypos.f0')-0.6
  dim=size(xp)
;
  fzread,dd,'../../../2017.07.09/calib/dd1/mihi_sp.dd1.sum.f0',h
  dd2=double(dd)/400.0
;
  fzread,dd,'../dd1/mihi_sp.dd1.sum.f0',h
  dd=double(dd)/400.0
;
 
files=['../ff1/mihi_sp.ff1.sum.f0','../../../2017.07.09/calib/ff1/mihi_sp.ff1.sum.f0','../ff2/mihi_sp.ff2.sum.f0','../../../2017.07.09/calib/ff2/mihi_sp.ff2.sum.f0',$
       
'../ff4/mihi_sp.ff4.sum.f0','../../../2017.07.09/calib/ff3/mihi_sp.ff3.sum.f0','../ff4/mihi_sp.ff4.sum.f0','../../../2017.07.09/calib/ff4/mihi_sp.ff4.sum.f0',$
       
'../ff5/mihi_sp.ff5.sum.f0','../../../2017.07.09/calib/ff5/mihi_sp.ff5.sum.f0','../ff6/mihi_sp.ff6.sum.f0',$
         '../../../2017.07.08/calib/ff3/mihi_sp.ff3.sum.f0']
 
files=['../ff1/mihi_sp.ff1.sum.f0','../ff2/mihi_sp.ff2.sum.f0','../ff4/mihi_sp.ff4.sum.f0','../ff5/mihi_sp.ff5.sum.f0','../ff6/mihi_sp.ff6.sum.f0',$
       
'../../../2017.07.09/calib/ff1/mihi_sp.ff1.sum.f0','../../../2017.07.09/calib/ff2/mihi_sp.ff2.sum.f0',$
       
'../../../2017.07.09/calib/ff3/mihi_sp.ff3.sum.f0','../../../2017.07.09/calib/ff4/mihi_sp.ff4.sum.f0',$
         '../../../2017.07.09/calib/ff5/mihi_sp.ff5.sum.f0']
;
files=['../ff1/mihi_sp.ff1.sum.f0','../ff2/mihi_sp.ff2.sum.f0','../ff4/mihi_sp.ff4.sum.f0','../ff5/mihi_sp.ff5.sum.f0','../ff6/mihi_sp.ff6.sum.f0']
;
  ff=dblarr(dim[2],dim[1],dim[3],3*n_elements(files))
  msk=dblarr(dim[2],dim[1],dim[3])
  for f=0,n_elements(files)-1 do begin
    fzread,flt,files[f],h
    dx=dd
    if(f>5) then dx=dd2
    ffx=map(double(flt)/10000.0-dx,xp,yp)
    help,ffx,ff
    med=median(ffx)
    for i=0,dim[3]-2 do begin
      xx=ffx[*,*,i+1]/median(ffx[*,*,i+1])
      mm=xx*0.0+1.0
      idx=where((xx gt 2.0)or(xx lt 0.5))
      mm[idx]=0.0
      msk[*,*,i]+=mm
      xx[idx]=0.0
      ff[*,*,i,3*f]=xx
    end
    for i=0,dim[3]-1 do begin
      xx=ffx[*,*,i]/median(ffx[*,*,i])
      mm=xx*0.0+1.0
      idx=where((xx gt 2.0)or(xx lt 0.5))
      mm[idx]=0.0
      msk[*,*,i]+=mm
      xx[idx]=0.0
      ff[*,*,i,3*f+1]=xx
    end
    for i=1,dim[3]-1 do begin
      xx=ffx[*,*,i-1]/median(ffx[*,*,i-1])
      mm=xx*0.0+1.0
      idx=where((xx gt 2.0)or(xx lt 0.5))
      mm[idx]=0.0
      msk[*,*,i]+=mm
      xx[idx]=0.0
      ff[*,*,i,3*f+2]=xx
    end
  end
  msk_idx=where(msk eq 3*n_elements(files))
;
  print,'m=',string(n_elements(msk_idx),format='(I8)'),'( out
of'+string(n_elements(ff[*,*,*,0]),format='(I8)')+')','x',n_elements(files)
;  m=transpose(reform(ff,[n_elements(ff[*,*,*,0]),n_elements(files)]))
  m=dblarr(3*n_elements(files),n_elements(msk_idx))
  for f=0,3*n_elements(files)-1 do begin
    ffn=ff[*,*,*,f]
    m[f,*]=ffn[msk_idx]
  end
  stats,m
  help,m
  svdc,m,w,u,v
  usz=size(u)
  for f=0,usz[1]-1 do print,sqrt(total(u[f,*]*u[f,*]))
;
;
;
  if(0 gt 1) then begin
    names=file_search('/SST/2017.07.10/obs07/mihi_sp/','*.f0')
    a=dd*0.0
    nn=0
    for i=(108900L)/3L,(109799L)/3L do begin
      a+=map(double(f0(names[i]))-dd,xp,yp)
      nn+=1
    end
    a/=nn
    fzwrite,a,'myavg.f0','bla'
  end else begin
    fzread,a,'myavg.f0',h
  end
;
  b=a
  for i=0,dim[3]-1 do a[*,*,i]=a[*,*,i]/median(a[*,*,i])

  cmx=usz[1]-1
  cfs=dblarr(cmx+1)
  for i=0,cmx do cfs[i]=total(a[msk_idx]*u[i,*])
  print,cfs
;
  ffc=a*0.0
  for i=0,cmx do ffc[msk_idx]+=cfs[i]*u[i,*]
;
  for i=0,dim[3]-1 do ffc[*,*,i]=median(ffc[*,*,i])/ffc[*,*,i]
  ffc[where((ffc gt 2.)or(ffc lt 0.5))]=0
  stats,ffc

  tvscl,histo_opt(b[*,*,100])
  tvscl,histo_opt(ffc[*,*,100]),128,0

  a=b*ffc
  tvscl,histo_opt(rebin((a)[*,*,100],228,256,/sample),0.01),256,0
stop
  names=file_search('/SST/2017.07.10/obs07/mihi_sp/','*.f0')
;  for k=56700L,183599L,900L do begin
  for k=133200L,183599L,900L do begin
    x=c*0.0
    nn=0
    for j=k/3L,(k+899L)/3L do begin
      a=map(double(f0(names[j]))-dd,xp,yp)
      b=a
      for i=0,dim[3]-1 do a[*,*,i]/=median(a[*,*,i])
;
      for i=0,cmx do cfs[i]=total(a[msk_idx]*u[i,*])
      print,cfs
      ffc=a*0.0
      for i=0,cmx do ffc[msk_idx]+=cfs[i]*u[i,*]
      for i=0,dim[3]-1 do ffc[*,*,i]=1.0/ffc[*,*,i]
      ffc[where((ffc gt 2.)or(ffc lt 0.5))]=0
      for i=0,dim[3]-1 do a[*,*,i]=a[*,*,i]*ffc[*,*,i]+median(b[*,*,i])
; 
      tvscl,histo_opt(rebin(a[*,*,140],228,256,/sample),0.01),512,0
      x+=a
      nn+=1

      bla=strsplit(names[j],'/',/extract)
    
fzwrite,float(a),'/arch/disk0/reduc/2017.07.10/obs07/ff4/'+bla[n_elements(bla)-1],'header'
      tvscl,histo_opt(rebin((x/nn)[*,*,140],228,256,/sample),0.01),768,0
    end
    x/=nn
  
fzwrite,float(x),'/arch/disk0/reduc/2017.07.10/obs07/ff4_avg/average.'+string(k,format='(I06)')+'..'+string(k+899L,format='(I06)')+'.f0','header'
  end
 
  stop
 
end


;; stop
;;
;; d=c
;; for i=0,dim[3]-1 do d[*,*,i]=d[*,*,i]/median(d[*,*,i])-1.0
;; cfs2=dblarr(cmx+1)
;; for i=0,cmx do cfs2[i]=total(d*u[i,*])/sqrt(total(u[i,*]*u[i,*]))
;; print,cfs2
;
;; ffc=d*0.0+1.0
;; for i=0,cmx do ffc+=cfs2[i]*u[i,*]
;
;; for i=0,dim[3]-1 do ffc[*,*,i]=1.0/ffc[*,*,i]
;  ffc[where((ffc gt 2.)or(ffc lt 0.5))]=0
;; stats,ffc
;; x=c*ffc

;; tvscl,histo_opt(rebin((x)[*,*,140],228,256,/sample),0.01),512,0

;; y=x
;; for i=0,dim[3]-1 do y[*,*,i]=y[*,*,i]/median(y[*,*,i])-1.0
;; cfs3=dblarr(cmx+1)
;; for i=0,cmx do cfs3[i]=total(y*u[i,*])/sqrt(total(u[i,*]*u[i,*]))
;; print,cfs3
;; ffc=d*0.0+1.0
;; for i=0,cmx do ffc+=cfs2[i]*u[i,*]
;
;; for i=0,dim[3]-1 do ffc[*,*,i]=1.0/ffc[*,*,i]
;  ffc[where((ffc gt 2.)or(ffc lt 0.5))]=0
;; stats,ffc
;; y=x*ffc

;; tvscl,histo_opt(rebin((y)[*,*,140],228,256,/sample),0.01),768,0

;; stop
;  splitlist=['xad','xba']
;
splitlist=['xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal','xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay','xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg']
;
flatlist=['../ff1/mihi_sp.ff1','../ff2/mihi_sp.ff2','../ff3/mihi_sp.ff3','../ff4/mihi_sp.ff4','../ff5/mihi_sp.ff5','../ff6/mihi_sp.ff6',$
;          
'../../../2017.07.09/calib/ff1/mihi_sp.ff1','../../../2017.07.09/calib/ff2/mihi_sp.ff2','../../../2017.07.09/calib/ff3/mihi_sp.ff3','../../../2017.07.09/calib/ff4/mihi_sp.ff4','../../../2017.07.09/calib/ff5/mihi_sp.ff5']
;          
'../../../2017.07.08/calib/ff1/mihi_sp.ff1','../../../2017.07.08/calib/ff2/mihi_sp.ff2','../../../2017.07.08/calib/ff3/mihi_sp.ff3','../../../2017.07.08/calib/ff4/mihi_sp.ff4','../../../2017.07.08/calib/ff5/mihi_sp.ff5']
;
flatlist=['../ff2/mihi_sp.ff2','../ff3/mihi_sp.ff3','../ff4/mihi_sp.ff4','../ff5/mihi_sp.ff5']
;  files=strarr(n_elements(flatlist)*n_elements(splitlist))
;  for f=0,n_elements(flatlist)-1 do for l=0,n_elements(splitlist)-1 do
files[f*n_elements(splitlist)+l]=flatlist[f]+'.'+splitlist[l]+'.sum.f0'