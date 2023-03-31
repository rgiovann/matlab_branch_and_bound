      % *****  MESTRANDO   : GIOVANNI LEOPOLDO ROZZA *************
      %        DISCIPLINA  : PROGAMAÇÃO  LINEAR INTEIRA E OTIMIZAÇÃO DE
      %                      REDES
      %                PPGEP - TURMA 2012
      % **********************************************************
      %  PROBLEMA ROTEIRIZACAO (CAIXEIRO VIAJANTE)
      %  SOLUÇÃO ATRAVÉS ALGORITMO DE DESIGNAÇÃO + BRANCH AND BOUND
      %  USANDO SOLVER CEPLEX CHAMADO ATRAVES DE API PELO MATLAB
      % **********************************************************
      
      global MAXIMO_NUMERO_RAMOS;
      global Aeq
      global beq;
      global Aineq
      global bineq;
      global matrizCustos_Atual;
      global na;
      global nb;
      global VALOR_INFINITO;
      global idx_NO_INICIAL;
      global idx_NIVEL_FOLHA;
      global ptr_FOLHA_MAE;
      global idx_VALOR_ATUAL_Z;
      
      clc; % clear command window
      
      MAXIMO_NUMERO_RAMOS = 1000;
      VALOR_INFINITO = 9999;
      MENOR_VALOR_ATUAL = 9999;
      CIRCUITO_MENOR_CUSTO=zeros(1,(na+1));
      matrizRestricoes=zeros(na*nb,2);
      % preciso armazenar todos os subcircuitos para imprimir!!
      matriz_Todos_Circuitos = zeros(1000,(na+2)); % ultimo campo é a folha mae
      ptr_TODOS_CIRCUITOS=0;
      
      % set CPLEX path
      addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Preview124\cplex\matlab\x64_win64')
      
      % obtendo a matriz de dados da Struct via arquivo Excel
      structDadosEntrada = uiimport('-file');
      arrayNomeDoCampoStruct = fieldnames(structDadosEntrada);
      matrizCustos_Inicial = structDadosEntrada.(arrayNomeDoCampoStruct{1,1});  
      
      [na,nb] = size(matrizCustos_Inicial);

      
      %************************************************************
      % armazena as folhas da árvore
      %************************************************************
      
      % armazeno os ramos em uma matriz de 1000 linhas, minha capacidade
      % portanto de armazenar o ramo é determinada pelo tamanho dessa
      % matriz. As primeiras n colunas são o tamanho (fixo em na+1) maximo 
      % do circuito.
      % A nésima+1 coluna é o nó i inicial que se conecta a todos os outros nós
      % para definir as restriçoes Cij=inf, 
      % o nésimo+2 é o nível da folha e
      % o nésimo+3 é ponteiro com o indice da folha-mãe
      % o nésimo+4 é o valor de Z (se igual a zero, nó ainda tem
      % subdivisoes para niveis inferiores), se = -1 entao todas os ramos
      % abaixo desse nível já tiveram seu bound (circuito com nr nós=na) encontrado.
      %
      %     n tam max do circuito
      %  |<------------------------->| n+1 | n+2 | n+3 | n+4 |
      %
      matrizFolhas = zeros(MAXIMO_NUMERO_RAMOS,(na+1) + 4 ) ;
      
      idx_NO_INICIAL     = (na+1) + 1;
      idx_NIVEL_FOLHA    = (na+1) + 2;
      ptr_FOLHA_MAE      = (na+1) + 3;
      idx_VALOR_ATUAL_Z  = (na+1) + 4;
             
      matrizCustos_Atual = matrizCustos_Inicial;    
      
      % **********************************************************
      % montando as matrizes de restrições (NUNCA MUDA)
      % **********************************************************
      
      Aineq=[];  % inequações nulas
      bineq=[];
            
      % matriz de restrições vai ter 2*i linhas e i*j colunas
      
      Aeq = zeros(2*na,na*nb);
      
      % inicilizo linha 1 e linha 1+na
      for i=1:na
          
          % somatoria de Akj = 1 
          Aeq(1:1,i:i)=1;  
          
          % somatorio de Aik = 1
          Aeq(1+na:1+na,1+(na*(i-1)):1+(na*(i-1)))=1;
          
      end
            
      % inicializo demais linhas fazendo shift "binário" de na posições nas
      % primeiras na restrições (linha constante) ou de 1 posição
      % nas na+1 posições (coluna constante)
      
      for i=2:na
          Aeq(i:i,:) = circshift(Aeq(1:1,:),[0,(i-1)*na]);
          Aeq(i+na:i+na,:) = circshift(Aeq(1+na:1+na,:),[0,(i-1)]);
      end
        
      % inicilizo o vetor de capacidades (matriz coluna) com valor 1 
      % (NUNCA MUDA)
      beq=ones(2*na,1);
      
      % **********************************************************
      % fim da montagem das matrizes de capacidade e restricoes
      % **********************************************************
      
      % **************************************
      % Chamo CPLEX solver atraves de funcao
      % *************************************  
      
      [matrizSolucao, nValorAtualZ] = ChamaSolverCPLEX();
           
      % **************************************
      % Extraio os subcircuitos
      % ************************************* 
      
      [ matrizCircuitosAtual ]  = MontaCircuitosdaSolucao(matrizSolucao);
            
      %*******************************************************      
      % descobre o menor circuito 
      %*******************************************************
      
      [menorCircuito] = AchaMenorCircuito(matrizCircuitosAtual);
      
      %*******************************************************      
      % atualiza as folhas do nivel 001
      %*******************************************************
      NIVEL_DO_RAMO_ATUAL=1;
      % armazena o menor circuito na matriz de folhas
      % armazena os nos de partida para montar as restrições
      
      lTerminouArvore=false;

      % se tem subcircuito Z = -1
      if find( menorCircuito(1:1,:),1,'last')  < na+1
         for i=1:(find( menorCircuito(1:1,:),1,'last')-1)
            matrizFolhas(i:i,1:na+1) = menorCircuito(1:1,1:na+1);
            matrizFolhas(i:i,idx_NO_INICIAL :idx_NO_INICIAL) = matrizFolhas(1:1,i:i);
            matrizFolhas(i:i,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA)= 1 ;
            matrizFolhas(i:i,ptr_FOLHA_MAE:ptr_FOLHA_MAE) = 0;  
            matrizFolhas(i:i,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z) = nValorAtualZ;         
         end
      
         ptr_MATRIZ_FOLHA_LAST = i;
         ptr_MATRIZ_RESTRICOES=0;
         % indice para viajar sobre o ponteiro
         % processa todas as folhas do ramo atual
      
         % armazena o ponteiro da nivel da folha corrente (working) onde serão
         % calculados todas as soluçoes
         ptr_NIVEL_FOLHA_WORK=ptr_MATRIZ_FOLHA_LAST;      
      else % nem executa proximas iteracoes, solucao otima na 1a designacao
          disp('**** Solução Ótima na primeira iteração ****');
          disp(menorCircuito);
          strLinha='Valor Z* =..';
          strLinha=strcat(strLinha,sprintf('%3.3d',nValorAtualZ));
          disp(strLinha);
          return
      end
      i=1;
      % armazena circuitos para imprimir depois
      [imCA,jmCA] =size(matrizCircuitosAtual);
      while i <= imCA
          if sum(matrizCircuitosAtual(i:i,1:(na+1))) >0
             ptr_TODOS_CIRCUITOS = ptr_TODOS_CIRCUITOS + 1;
             matriz_Todos_Circuitos(i:i,1:(na+1)) = matrizCircuitosAtual(i:i,1:(na+1));
             matriz_Todos_Circuitos(i:i,(na+2):(na+2)) = 0;
          end
          i=i+1;
      end
      while ~lTerminouArvore
            
            %*********************************                
            % monta restriçoes anteriores
            %*********************************
            % lê o ponteiro da folha-mae
            ptr_Folha = matrizFolhas(ptr_NIVEL_FOLHA_WORK:ptr_NIVEL_FOLHA_WORK,ptr_FOLHA_MAE:ptr_FOLHA_MAE);
            while true
               if ptr_Folha ~= 0
                  nPontoInicial = matrizFolhas(ptr_Folha:ptr_Folha,idx_NO_INICIAL:idx_NO_INICIAL);
                  menorCircuito = matrizFolhas(ptr_Folha:ptr_Folha,1:(na+1));
                  k=find( menorCircuito(1:1,:),1,'last')-1;  % nao leva em consideracao ultimo nó que é igual ao primeiro
                  for m=1:k
                     if matrizFolhas(ptr_Folha:ptr_Folha,m:m) ~= nPontoInicial
                        ptr_MATRIZ_RESTRICOES = ptr_MATRIZ_RESTRICOES +1;
                        matrizRestricoes(ptr_MATRIZ_RESTRICOES:ptr_MATRIZ_RESTRICOES,1:1) = nPontoInicial;
                        matrizRestricoes(ptr_MATRIZ_RESTRICOES:ptr_MATRIZ_RESTRICOES,2:2) = matrizFolhas(ptr_Folha:ptr_Folha,m:m);
                     end
                  end 
               
                  % busca proxima folha-mae, se for zero, chegou na
                  % raiz
                  ptr_Folha = matrizFolhas(ptr_Folha:ptr_Folha,ptr_FOLHA_MAE:ptr_FOLHA_MAE);             
               else
                  % folha-mae é zero, sai do loop while..end
                  break;
               end
               % 
            end
            
            %*********************************                
            % monta restriçoes da folha corrente
            %*********************************
            nPontoInicial = matrizFolhas(ptr_NIVEL_FOLHA_WORK:ptr_NIVEL_FOLHA_WORK,idx_NO_INICIAL:idx_NO_INICIAL);
            menorCircuito = matrizFolhas(ptr_NIVEL_FOLHA_WORK:ptr_NIVEL_FOLHA_WORK,1:(na+1));
            k=find( menorCircuito(1:1,:),1,'last')-1;  % nao leva em consideracao ultimo nó que é igual ao primeiro
            for m=1:k
               if menorCircuito(1:1,m:m) ~= nPontoInicial
                  ptr_MATRIZ_RESTRICOES = ptr_MATRIZ_RESTRICOES +1;
                  matrizRestricoes(ptr_MATRIZ_RESTRICOES:ptr_MATRIZ_RESTRICOES,1:1) = nPontoInicial;
                  matrizRestricoes(ptr_MATRIZ_RESTRICOES:ptr_MATRIZ_RESTRICOES,2:2) = matrizFolhas(ptr_NIVEL_FOLHA_WORK:ptr_NIVEL_FOLHA_WORK,m:m);
               end               
            end
            
            %*********************************
            % monta as restricoes Cij=infinito 
            %*********************************
            % le a matriz de custos original
            matrizCustos_Atual = matrizCustos_Inicial;
            for m=1:ptr_MATRIZ_RESTRICOES
               i=matrizRestricoes(m:m,1:1);
               j=matrizRestricoes(m:m,2:2);
               matrizCustos_Atual( i:i,j:j) = VALOR_INFINITO;         
            end
            
            %resolve problema designação.
             
            % **************************************
            % Chamo CPLEX solver atraves de funcao
            % *************************************  
      
            [matrizSolucao, nValorAtualZ] = ChamaSolverCPLEX();
           
            % **************************************
            % Extraio os subcircuitos
            % ************************************* 
      
            [ matrizCircuitosAtual ] = MontaCircuitosdaSolucao(matrizSolucao);
            
            %*******************************************************      
            % descobre o menor circuito 
            %*******************************************************
      
            [menorCircuito] = AchaMenorCircuito(matrizCircuitosAtual);
            
            % armazena circuitos para imprimir depois
            i=1;
            [imCA,jmCA] =size(matrizCircuitosAtual);
            while i<=imCA
               if sum(matrizCircuitosAtual(i:i,1:(na+1))) >0
                  ptr_TODOS_CIRCUITOS = ptr_TODOS_CIRCUITOS + 1;
                  matriz_Todos_Circuitos(ptr_TODOS_CIRCUITOS:ptr_TODOS_CIRCUITOS,1:(na+1)) = matrizCircuitosAtual(i:i,1:(na+1));
                  matriz_Todos_Circuitos(ptr_TODOS_CIRCUITOS:ptr_TODOS_CIRCUITOS,(na+2):(na+2)) = ptr_NIVEL_FOLHA_WORK;
               end
               i=i+1;
            end
                        
            %*******************************************************
            % atualiza a matriz de folhas (nova folha criada)            
            % armazena o menor circuito na matriz de folhas
            % armazena os nos de partida para montar as restrições caso 
            % tenha subcircuito. Se tem subcircuito Z = -1
            %*******************************************************
            if find( menorCircuito(1:1,:),1,'last')  < na+1
                if nValorAtualZ < MENOR_VALOR_ATUAL % menor que o valor otimo, processa restricoes
                   for i=1:(find( menorCircuito(1:1,:),1,'last')-1)
                      ptr_MATRIZ_FOLHA_LAST = ptr_MATRIZ_FOLHA_LAST + 1;              
                      % atualiza ponteiro da matriz de folhas
                      matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,1:na+1) = menorCircuito(1:1,1:na+1);
                      matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NO_INICIAL :idx_NO_INICIAL) = menorCircuito(1:1,i:i);
                      % atualiza para novo nivel 
                      matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA)= NIVEL_DO_RAMO_ATUAL + 1 ;
                      % armazena indice da folha-mae
                      matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,ptr_FOLHA_MAE:ptr_FOLHA_MAE) = ptr_NIVEL_FOLHA_WORK;      
                      % nao tem valor de Z pois é subcircuito
                      matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z) = nValorAtualZ;         
                   end
                else % maior que o valor otimo, nao faz sentido processar restricoes (mata a folha).
                   ptr_MATRIZ_FOLHA_LAST = ptr_MATRIZ_FOLHA_LAST + 1;              
                   matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,1:na+1) = menorCircuito(1:1,1:na+1);
                   % nao tem no inicial, pois é folha com Z maior que o otimo atual(folha morre aqui) 
                   matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NO_INICIAL :idx_NO_INICIAL) = -1;
                   % nao tem novo nivel pois é folha com Z maior que o otimo atual (folha morre aqui) 
                   matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA)= -1 ;
                   % armazena indice da folha-mae
                   matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,ptr_FOLHA_MAE:ptr_FOLHA_MAE) = ptr_NIVEL_FOLHA_WORK;      
                   % armazena valor de Z
                   matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z) = nValorAtualZ;    
                end
            else % achou circuito completo!
                 ptr_MATRIZ_FOLHA_LAST = ptr_MATRIZ_FOLHA_LAST + 1;              
                 matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,1:na+1) = menorCircuito(1:1,1:na+1);
                 % nao tem no inicial, é circuito completo
                 matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NO_INICIAL :idx_NO_INICIAL) = -1;
                 % nao tem novo nivel pois é circuito completo (nivel da folha
                 % mae)
                 matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA)= -1 ;
                 % armazena indice da folha-mae
                 matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,ptr_FOLHA_MAE:ptr_FOLHA_MAE) = ptr_NIVEL_FOLHA_WORK;      
                 % armazena valor de Z
                 matrizFolhas(ptr_MATRIZ_FOLHA_LAST:ptr_MATRIZ_FOLHA_LAST,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z) = nValorAtualZ;    
                 % verifica se é a melhor solução até agora
                 if nValorAtualZ < MENOR_VALOR_ATUAL
                     MENOR_VALOR_ATUAL    = nValorAtualZ; 
                     CIRCUITO_MENOR_CUSTO = menorCircuito; 
                 end
            end                   
            
            %*****************************************
            % passo para a folha anterior do nível atual
            %*****************************************
            
            % reinicializa matriz restricoes e ponteiro de restricoes
            matrizRestricoes=zeros(na*nb,2);
            ptr_MATRIZ_RESTRICOES=0;
            
            ptr_NIVEL_FOLHA_WORK = ptr_NIVEL_FOLHA_WORK - 1;
            
            % verifica se tem novos niveis ou ja processou tudo, assume que
            % sim
            lTerminouArvore=true;  
            
            if ptr_NIVEL_FOLHA_WORK == 0  % cheguei no fim do 1o nivel
                % atualiza ponteiro para o  ultimo elemento do proximo
                % nivel que sera agora corrente
                
                % 1o nivel criou novos niveis ?
                i=ptr_MATRIZ_FOLHA_LAST;
                while i > 0
                   k= matrizFolhas(i:i,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA);
                   
                   % existe novo nivel ?, senao termina programa
                   if k == NIVEL_DO_RAMO_ATUAL +1
                       lTerminouArvore=false; 
                       NIVEL_DO_RAMO_ATUAL = NIVEL_DO_RAMO_ATUAL +1;
                       break;
                   end
                   i=i-1;
                end
                
                if ~lTerminouArvore
                    
                    % ainda tem mais um nivel a processar, posiciona o
                    % ponteiro na primeira folha que tenha subcircuitos
                   i=ptr_MATRIZ_FOLHA_LAST;
                   while true
                       k= matrizFolhas(i:i,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA);
                       if k ~= -1
                          break;
                       end
                       i=i-1;
                   end
                   ptr_NIVEL_FOLHA_WORK = i;
                end
            else
                
                 % não é + o fim do 1o nivel,entao avalia se existem outroas folhas com o
                 % mesmo nível para serem processados, filtra os que tem
                 % nivel = -1 (nao tem subcircuito)
                 
                 i=ptr_NIVEL_FOLHA_WORK;
                
                 while true
                    k= matrizFolhas(i:i,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA);
                    if k == -1
                       % circuito completo, continua busca
                       i=i-1;
                    else
                       if k == NIVEL_DO_RAMO_ATUAL
                          % ainda tem folha para processar, sai do loop
                          lTerminouArvore=false; 
                          ptr_NIVEL_FOLHA_WORK = i;
                          break
                       else
                           % entrou no nivel anterior 
                           % ja processado, sai do loop
                           break;
                       end
                    end
                 end
                 
                 
                 % se não achou mais linhas com mesmo nivel, já processou
                 % tudo, logo devo verificar se novos niveis foram criados
                 % ou as folhas processadas nao geraram mais subcircuitos
                 
                 if lTerminouArvore
                     i = ptr_MATRIZ_FOLHA_LAST;
                     while i > ptr_NIVEL_FOLHA_WORK
                        k= matrizFolhas(i:i,idx_NIVEL_FOLHA:idx_NIVEL_FOLHA);  
                        if k == NIVEL_DO_RAMO_ATUAL +1
                           lTerminouArvore = false; 
                           % tem mais um novo nivel para processar, atualiza
                           % ponteiro de nivel de trabalho
                           ptr_NIVEL_FOLHA_WORK = i ;
                           NIVEL_DO_RAMO_ATUAL  = NIVEL_DO_RAMO_ATUAL + 1;
                           break;
                        end
                        i=i-1;
                     end
                 end
            end                       
         
      end
        
      
      % *****************************
      % IMPRIME EVOLUCAO DO ALGORITMO
      % *****************************
      disp('*************************** ++ MATRIZ DE CUSTOS ORIGINAL ++ ***************************');
      disp(matrizCustos_Inicial);
      disp('**************************************************************************************');
      disp('.');
      disp('.');      
      disp('***************************** ++ SOLUÇÃO INICIAL ++ *****************************');
      j=ptr_TODOS_CIRCUITOS;
      str_subcircuitos='CIRCUITO(S) ->';
      while j > 0
         if matriz_Todos_Circuitos(j:j,(na+2):(na+2)) == 0  
            k=find( matriz_Todos_Circuitos(j:j,1:(na+1)),1,'last');
            str_subcircuitos=strcat(str_subcircuitos, '[');
            for m=1:k
               if m < k
                  str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',matriz_Todos_Circuitos(j:j,m:m)),'-');
               else
                  str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',matriz_Todos_Circuitos(j:j,m:m)));
               end
            end
            str_subcircuitos=strcat(str_subcircuitos, ']');
         end
         j=j-1;
      end
      disp(str_subcircuitos);
      strLinha='VALOR DE Z INICIAL ->';
      strLinha= strcat(strLinha,sprintf('%4.3d',matrizFolhas(1:1,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z)));
      disp(strLinha);
      disp('**************************************************************************************');
      disp('.');
      disp('.'); 
      disp('-----------------------------------------------------------------------------------------')
      disp('FL | MAE | CIRCUITO                      | RESTRIÇÕES=INF |  Z |');
      disp('-----------------------------------------------------------------------------------------')
      % imprime sub circuitos
      % procura a folha mae (lembre que estamos gerando novas folhas,
      % devemos imprimir as folhas ja criadas que estao sendo
      % processadas
      ptr_MTX_FOLHA=1;
      while sum(matrizFolhas(ptr_MTX_FOLHA:ptr_MTX_FOLHA,:)) > 0  
         strLinha='';
         str_subcircuitos='';
         k=matrizFolhas(ptr_MTX_FOLHA:ptr_MTX_FOLHA,ptr_FOLHA_MAE:ptr_FOLHA_MAE);
         if k~=0
            strLinha = strcat(strLinha,sprintf('%3.0d',ptr_MTX_FOLHA),'|',sprintf('%5.0d',k),'|');
         else
           strLinha = strcat(strLinha,sprintf('%3.0d',ptr_MTX_FOLHA),'| INIC ','|');
         end
         % subcircuitos com a mesma folha mae
         j=ptr_TODOS_CIRCUITOS;
         while j > 0
            if matriz_Todos_Circuitos(j:j,(na+2):(na+2)) == ptr_MTX_FOLHA  
               k=find( matriz_Todos_Circuitos(j:j,1:(na+1)),1,'last');
               str_subcircuitos=strcat(str_subcircuitos, '[');
               for m=1:k
                  if m < k
                     str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',matriz_Todos_Circuitos(j:j,m:m)),'-');
                  else
                     str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',matriz_Todos_Circuitos(j:j,m:m)));
                  end
               end
               str_subcircuitos=strcat(str_subcircuitos, ']');
            end
            j=j-1;
         end
      
         strLinha=strcat(strLinha,str_subcircuitos);
         % imprime restricao da folha corrente
         str_restricoes='';  
         nPontoInicial = matrizFolhas(ptr_MTX_FOLHA:ptr_MTX_FOLHA,idx_NO_INICIAL:idx_NO_INICIAL);
         menorCircuito = matrizFolhas(ptr_MTX_FOLHA:ptr_MTX_FOLHA,1:(na+1));
         k=find( menorCircuito(1:1,:),1,'last')-1;  % nao leva em consideracao ultimo nó que é igual ao primeiro
         for m=1:k
            if menorCircuito(1:1,m:m) ~= nPontoInicial
               str_restricoes = strcat(str_restricoes, 'C',sprintf('%1.0d',nPontoInicial),sprintf('%1.0d',menorCircuito(1:1,m:m)),'.' );
            end               
         end          
         strLinha=strcat(strLinha,' | ',str_restricoes);
         i=ptr_MATRIZ_FOLHA_LAST;
         lFim=true;
         while i > 0
             if matrizFolhas(i:i,ptr_FOLHA_MAE:ptr_FOLHA_MAE) == ptr_MTX_FOLHA
                strLinha=strcat(strLinha,' | ',sprintf('%4.3d',matrizFolhas(i:i,idx_VALOR_ATUAL_Z:idx_VALOR_ATUAL_Z)));
                lFim=false;
                break;
             end   
             i=i-1;
         end
         if ~lFim  % não é uma linha sem folha mae (circuito inteiro)
            disp(strLinha);
         end
         ptr_MTX_FOLHA=ptr_MTX_FOLHA+1;
      end
      disp('.');
      disp('.'); 
      strLinha='VALOR ÓTIMO DE Z ->';
      strLinha= strcat(strLinha,sprintf('%4.3d',MENOR_VALOR_ATUAL));
      disp(strLinha);
      str_subcircuitos='CIRCUITO DE MENOR CUSTO -> ';
      k=find( CIRCUITO_MENOR_CUSTO(1:1,1:(na+1)),1,'last');
      str_subcircuitos=strcat(str_subcircuitos, '[');
      for m=1:k
         if m < k
            str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',CIRCUITO_MENOR_CUSTO(1:1,m:m)),'-');
         else
               str_subcircuitos= strcat(str_subcircuitos,sprintf('%0.0d',CIRCUITO_MENOR_CUSTO(1:1,m:m)));
         end
      end
      str_subcircuitos=strcat(str_subcircuitos, ']');
      disp(str_subcircuitos);
      disp('.');
      disp('.');     
      disp('***** FIM PROGRAMA *****');

      % *****************************
      % FINAL IMPRESSÃO
      % *****************************
      
      
%
% vai achando as solucoes para todas as folhas do mesmo ramo, usa vetor de
% restricoes para guardar a memoria das restricoes a medida que desce os
% ramos, vai gerando novos nos para o proximo nível. O corte e quando não
% existir nenhum nivel criado pelo nivel anterior.
%
%
%
%
%
%
%

      
      
      
     
      
