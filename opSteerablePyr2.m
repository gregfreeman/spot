classdef opSteerablePyr2 < opSpot
   %OPSTEERABLEPYR   Sreerable Pyramid operator.
   %
   %   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties( SetAccess = private, GetAccess = public )
      levels      = 5;              % Number of levels
      orientations= 6;              % Number of orientations
      pind                          % Band Indicies
      masks                          % Band Masks
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
      function op = opSteerablePyr2(p,q,levels,orientations)
         
        n=prod([p q]);  
        

        [masks,pind] = opSFpyrInit([p,q],levels,orientations-1);
        m=sum(pind(:,1).*pind(:,2));

         
        op = op@opSpot('SteerablePyr2', m, n);
        op.signal_dims = [p, q];
        op.levels = levels;
        op.pind=pind;
        op.masks=masks;
        op.orientations=orientations;
         
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
            [y] = opSFpyrBuild(Xmat,op.levels,op.orientations,op.pind,op.masks);
         
         else % mode == 2
            y = opSFpyrRecon(x,op.levels,op.orientations,op.pind,op.masks);
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
