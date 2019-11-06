##-----obtain the initial values of parameters in different distributions----
init.para<-function(x,r,p,alpha1,alpha2,n){
r=r[r>0]
p=p[(p>0)&(p<1)]
alpha1=alpha1[alpha1>0]
alpha2=alpha2[alpha2>0]
n=n[n>0]
if(length(r)>0){
   r=r[1]  ##only use the first value if user inputs more than one applicable values
}else{
   r=max(x)  ##if there is no r input, use max(x) as initial value.
}
if(length(p)>0){
   p=p[1]  ##only use the first value if user inputs more than one applicable values
}else{
   p=sum(x>0)/length(x)  ##if there is no p input, use sum(x>0)/length(x) as initial value.
} 
if(length(n)>0){
   n=n[1]  ##only use the first value if user inputs more than one applicable values
}else{
   n=max(x)+1  ##if there is no p input, use max(x)+1 as initial value.
} 
temp1=mean(x)
temp2=stats::var(x)
if(length(alpha1)>0){
   alpha1=alpha1[1]  ##only use the first value if user inputs more than one applicable values
}else{
   alpha1=abs(temp1*(temp1*(1-temp1)/temp2-1))  ##the initial value derives from the moment estimate of beta distribution mean(mean(1-mean)/var-1).
} 
if(length(alpha2)>0){
   alpha2=alpha2[1]  ##only use the first value if user inputs more than one applicable values
}else{
   alpha2=abs((1-temp1)*(temp1*(1-temp1)/temp2-1))  ##the initial value derives from the moment estimate of beta distribution (1-mean)(mean(1-mean)/var-1).
} 
return(list(r=r,p=p,alpha1=alpha1,alpha2=alpha2,n=n))
}

