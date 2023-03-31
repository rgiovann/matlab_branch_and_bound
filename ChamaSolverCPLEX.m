function [ matrizDesignacao, nvalorZ ] = ChamaSolverCPLEX( ~ )


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

      %global MAXIMO_NUMERO_RAMOS;
      global Aeq
      global beq;
      global Aineq
      global bineq;
      global matrizCustos_Atual;
      global na;
      global nb;
            
      % formatando a matriz entrada em uma matriz coluna para CPLEX
      % função objetivo

      ObjFunction=zeros(na*nb,1);
      
      k=1;
      for i=1:na
          for j=1:nb
              ObjFunction(k:k,1:1) = matrizCustos_Atual(i:i,j:j);
              k=k+1;
          end
      end
      
      options = cplexoptimset;
      options.Diagnostics = 'off';
   
      [x, fval, exitflag, output] = cplexbilp (ObjFunction, Aineq, bineq, Aeq, beq);
      
      matrizSolucao = x;
      nvalorZ       = fval;
      
      %*********************************
      % monto a matriz de solução ixj
      % com as designações
      %*********************************
      
      matrizDesignacao = zeros(na,nb);
      k=1;
      for i=1:na
          for j=1:nb
              matrizDesignacao(i:i,j:j) = x(k:k,1:1);
              k=k+1;
          end
      end

end

