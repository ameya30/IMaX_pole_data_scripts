pro findstuff,numobse1,numcycl1

de =  where(numobse1 eq 300, kitna)
temp = []
foreach i, de DO BEGIN
  if numcycl1[i] ge 6 and numcycl1[i] le 28 then begin
    temp[i] = i
  endif
endforeach
end