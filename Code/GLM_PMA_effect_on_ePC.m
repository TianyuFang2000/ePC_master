%% Fig 3 plot: Edge-wise PMA effect on edge PC 
% ePC Development of FT
load('/data/users/tfang/eFC_entropy/eFC_entropy/Entropy_developing/Code/tFang_dHCPts_0314.mat')
indFT(144)=[];
FT_age=table2array(Info(indFT,"scan_age"));
FT_birth_age=table2array(Info(indFT,"birth_age"));
FT_gender=table2array(Info(indFT,"sex"));
Y=PC_FT';
for i=1:4005
    tbl1=table(FT_age,FT_birth_age,FT_gender,Y(:,i));
    glm1 = fitglm(tbl1);
    P(i)=table2array(glm1.Coefficients(2,4));    
    T_val(i)=table2array(glm1.Coefficients(2,3));   
    Beta(i)=table2array(glm1.Coefficients(2,1));  
    P_birth(i)=table2array(glm1.Coefficients(3,4));
end
% ePC Development of PT1 and PT2
PT1_age=table2array(Info(indPT1,"scan_age"));
PT1_birth_age=table2array(Info(indPT1,"birth_age"));
PT1_gender=table2array(Info(indPT1,"sex"));
Y1=PC_PT1_1';
PT2_age=table2array(Info(indPT2,"scan_age"));
PT2_birth_age=table2array(Info(indPT2,"birth_age"));
PT2_gender=table2array(Info(indPT2,"sex"));
Y2=PC_PT2_1';
for i=1:4005
    tbl1=table(PT1_age,PT1_birth_age,PT1_gender,Y1(:,i));
    glm1 = fitglm(tbl1);
    P1(i)=table2array(glm1.Coefficients(2,4));    
    T_val_PT1(i)=table2array(glm1.Coefficients(2,3));   
    tbl2=table(PT2_age,PT2_birth_age,PT2_gender,Y2(:,i));
    glm2 = fitglm(tbl2);
    P2(i)=table2array(glm2.Coefficients(2,4));    
    T_val_PT2(i)=table2array(glm2.Coefficients(2,3)); 
end