function [optValue,bestP,conv]=EHRJAYA(NP,Max_NFEs,LB,UB,Dim,fobj)

    addpath('../');
%     [LB ,UB,Dim ] = Parameter(pType);

%     index=pType;
    c=0.1;
    MaMi=(repmat((UB-LB),NP,1));
    MiLB=repmat((LB),NP,1);
    X=initialization(NP,Dim,UB,LB);
    NPg_current=NP;
    NPmin=3;
    Ra=0.3;
    uF1=0.5;
    uF2=0.5;

    for i=1:NP
        fitnessX(i)=fobj(X(i,:));
    end
    NFEs=NP;

    t=1;
    while NFEs<Max_NFEs
        Sf1=[];
        Sf2=[];
        n0=1;
        [fitnessBestX,~]=min(fitnessX);
        [fitnessWorstX,~]=max(fitnessX);
        [~,sortIndexX]=sort(fitnessX);
        Xbest=X(sortIndexX(1),:);
        Xworst=X(sortIndexX(end),:);
          for i=1:NPg_current
              R(sortIndexX(i))=i;
              P(sortIndexX(i)) = (R(sortIndexX(i))/NPg_current);
          end
          for i=1:NPg_current
            F1(i)=cauchyrnd(uF1,0.1);
            F2(i)=cauchyrnd(uF2,0.1);
            while (F1(i)<=0||F2(i)<=0)
                F1(i)=cauchyrnd(uF1,0.1);
                F2(i)=cauchyrnd(uF2,0.1);
            end
            if (F1(i)>1||F2(i)>1)
                F1(i)=1;
                F2(i)=1;
            end
        end

        for i=1:NPg_current
            if rand<0.5
        if fitnessBestX == 0
            W1=1;
        else
            W1=(mean(fitnessX)/fitnessBestX);
        end
        if fitnessWorstX == 0
            W2=1;
        else
            W2=(mean(fitnessX)/fitnessWorstX);
        end    
            rand1=rand;
            rand2=rand;
            
            for j=1:Dim 
                if P(i)<=1/3 && P(i)>0
                V(i,j) = X(i,j)+W1*rand1*(Xbest(j)-abs(X(i,j)))-W2*rand2*(Xworst(j)-abs(X(i,j)));
                elseif P(i)<=2/3 && P(i)>1/3
                XM=mean(X);
                V(i,j) = X(i,j)+W1*rand1*(Xbest(j)-abs(X(i,j)))-W2*rand2*(XM(j)-abs(X(i,j)));        
                else
                r1 = randi(NPg_current);
                  while r1==i
                      r1 = randi(NPg_current);
                  end
            r2 = randi(NPg_current);
                  while r2==r1 || r2==i
                      r2 = randi(NPg_current); 
                  end
                V(i,j) = X(i,j)+rand1*(Xbest(j)-abs(X(i,j)))+rand2*(X(r1,j)-X(r2,j));
                end
            end
            
            else
            r1 = randi(NPg_current);
                  while r1==i
                      r1 = randi(NPg_current);
                  end
            r2 = randi(NPg_current);
                  while r2==r1 || r2==i
                      r2 = randi(NPg_current); 
                  end
             for j=1:Dim 
               if P(i)<=1/3 && P(i)>0
                 V(i,j) = X(i,j)+F1(i).*(Xbest(j)-Xworst(j))+F2(i).*(Xbest(j)-X(r2,j));
               elseif P(i)<=2/3 && P(i)>1/3
                   XM=mean(X);
                   V(i,j) = X(i,j)+F1(i).*(Xbest(j)-Xworst(j))+F2(i).*(XM(j)-X(r2,j));
               else
                   V(i,j) = X(i,j)+F1(i).*(Xbest(j)-Xworst(j))+F2(i).*(X(r1,j)-X(r2,j));
               end
             end
            end
            for j=1:Dim
                if V(i,j)>UB || V(i,j)<LB
                    V(i,j)=LB+rand*(UB-LB);
                end
            end
            fitnessV(i)=fobj(V(i,:));
            NFEs=NFEs+1;
           
            if fitnessV(i)<fitnessX(i)
                Sf1(n0)=F1(i);
                Sf2(n0)=F2(i);
                n0=n0+1;
                X(i,:)=V(i,:);
                fitnessX(i)=fitnessV(i);
            end
        end
        
        [~,ab]=size(Sf1);
        if ab~=0
            newSf1=(sum(Sf1.^2))/(sum(Sf1));
            uF1=(1-c)*uF1+c.*newSf1;
            newSf2=(sum(Sf2.^2))/(sum(Sf2));
            uF2=(1-c)*uF2+c.*newSf2;
        end
        
        NPg_current=round(((NPmin-NP)/Max_NFEs)*NFEs+NP);
        numP=size(X,1);
        if NPg_current<numP
            [~,indexSortP2]=sort(fitnessX);
            X(indexSortP2(NPg_current+1:numP),:)=[]; 
            fitnessX(indexSortP2(NPg_current+1:numP))=[]; 
        end
        
        if rand<Ra 
            maxBoun=max(X);
            minBoun=min(X);
            for i=1:NPg_current
                k=2-2*NFEs/Max_NFEs; 
                GOP(i,:)=(((minBoun+maxBoun)/2)+(((minBoun+maxBoun))/2*k)-(X(i,:)/k)); 
                for j=1:Dim
                    if GOP(i,j)>maxBoun(j) || GOP(i,j)<minBoun(j)
                        GOP(i,j)=minBoun(j)+rand*(maxBoun(j)-minBoun(j));
                    end
                end
                fitnessGOP(i)=fobj(GOP(i,:));
                NFEs=NFEs+1;
            end
            fitnessNew=[fitnessX,fitnessGOP];
            NewP=[X;GOP];
            [sortFitnessNew,sortIndex]=sort(fitnessNew);
            X=NewP(sortIndex(1:NPg_current),:);
            fitnessX=sortFitnessNew(1:NPg_current);
        end
        [fitnessBestX,recordIndex]=min(fitnessX);
        conv(t)=fitnessBestX;
        t=t+1;
    end
   
    endNFEs=NFEs;
    optValue=fitnessBestX;
    bestP=X(recordIndex,:);

end
