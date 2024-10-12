function [CL,CD,Cl,Cm,Cn] = avl(avlplane,operatingpoint)
%% Make a plane AVL file
avlplane.writeAVL('tempplane');

%% Make a op-point RUN file
%operatingpoint.writeRUN('tempplane');

%% Make a runfile
fh=fopen('temprunfile.txt','w');
fprintf(fh,'oper\n');
fprintf(fh,'a\n');
fprintf(fh,'a %f\n',operatingpoint.alpha);
fprintf(fh,'x\n');
fprintf(fh,'\n');
fprintf(fh,'quit');
fclose(fh);
if ispc
    command=sprintf('avl tempplane.avl < temprunfile.txt > tempoutput.txt');
elseif ismac
    command=sprintf('%s/avl3.35 %s/tempplane.avl < %s/temprunfile.txt > %s/tempoutput.txt',pwd,pwd,pwd,pwd);
end
system(command);

%% Analyze Runfile
rawtext=importdata('tempoutput.txt','',10000);
rawtext=string(rawtext);
for i=1:length(rawtext)
    if contains(rawtext(i),'Run case:')
        break;
    end
end
if i==length(rawtext)
    error()
end
interestingstuffstart=i+1;
interestingstuffend=interestingstuffstart+10;

data=rawtext(interestingstuffstart:interestingstuffend);

if contains(data(7),'CLtot')
    CL=double(extractAfter(data(7),'='));
else
    error()
end
if contains(data(8),'CDtot')
    CD=double(extractAfter(data(8),'='));
else
    error()
end
if contains(data(4),"Cltot")
    Cl=double(extractBetween(data(4),"Cltot =","Cl'tot"));
else
    error()
end
if contains(data(5),"Cmtot")
    Cm=double(extractAfter(data(5),"Cmtot ="));
else
    error()
end
if contains(data(6),"Cntot")
    Cn=double(extractBetween(data(6),"Cntot =","Cn'tot"));
else
    error()
end


end

