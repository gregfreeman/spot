classdef opSteerablePyr < opSpot
   %OPSTEERABLEPYR   Sreerable Pyramid operator.
   %
   %   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties( SetAccess = private, GetAccess = public )
      levels      = 5;              % Number of levels
      orientations= 6;              % Number of orientations
      pind                          % Indicies
      signal_dims                  % Dimensions of the signal domain
      funHandle                    % Multiplication function
      funHandle2                   % Divide function
   end % Properties
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opWavelet. Constructor.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opSteerablePyr(p,q,levels,orientations)
         
        n=prod([p q]);  
        

        [pyr,pind] = buildSFpyr(zeros(p,q),levels,orientations-1);
        m=numel(pyr);
         
         op = op@opSpot('SteerablePyr', m, n);
         op.signal_dims = [p, q];
         op.levels = levels;
         op.pind=pind;
       
         % Initialize function handle
         op.funHandle = @multiply_intrnl;
         op.funHandle2 = @divide_intrnl;
         
      end % function opWavelet
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Divide
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = mldivide(op,x)
         y = op.funHandle2(op,x);
      end % function multiply
      
   end % methods - public
   
   methods( Access = private )
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % matvec.  Application of Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply_intrnl(op,x,mode)
         p = op.signal_dims(1);
         q = op.signal_dims(2);
         if issparse(x), x = full(x); end

         if mode == 1
            Xmat = reshape(x,p,q);
            [y,~] = buildSFpyr(Xmat,op.levels,op.orientations-1);
         
         else % mode == 2
            y = reconSFpyr(x, op.pind);
            y=y(:);
         end
      end % function matvec
      

      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % divide_intrnl.  Application of redundant Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = divide_intrnl(op,x)       
          
         if issparse(x), x = full(x); end
         
         y = reconSFpyr(x, op.pind, op.levels, op.orientations-1);
         
      end % function divide      
         
   end % methods - private
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Multiply
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply(op,x,mode)
         y = op.funHandle(op,x,mode);
      end % function multiply
      
   end % methods - protected
   
end % classdef
