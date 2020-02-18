function h=errorshade(x,y,le,ue,c,error_est_type)
%%function h=errorshade(x,y,ue,le,c,error_est_type)
%error_est_type is a string: 'scalar' or 'vector'
%this function can return a handle for manipulation of graphics properties 
%or can be called sans handle and will plot shaded error bars 
%x = vector x data, y = vector y data, ue and le = upper/lower y error values
%respectively, c = color (input: 'c')
%if want symmetric error bars, enter in same vector or number for ue and le 
%(ie can enter vector so error bar is diff at each value or can enter in a 
%number, eg 95% conf int low and up bound (same number, one does not have 
%to be negative) and you'll get shaded in error bars for cont data)
%RS 

x=x(:);y=y(:);le=le(:);ue=ue(:);

if strcmp(error_est_type,'scalar')==1
thebottom=y-le;
%thetop=y+e;
thetop=y+ue;
else
    thebottom=le;
    thetop=ue;
end;

theregiony=[thebottom(1:end);thetop(end:-1:1)];
theregionx=[x(1:end);x(end:-1:1)];
h=patch(theregionx,theregiony,c,'Edgecolor',c,'FaceAlpha',0.2,'EdgeAlpha',0.15); %.2 = norm for big graph


%e.g.
%figure(1)
%errorshade(x,y,ue,le,c)
%hold on
%plot(x,y,'r','LineWidth',1.2)  



%SEM and Conf Int in matlab: http://www.matlab-cookbook.com/recipes/0100_Statistics/010_sem.html