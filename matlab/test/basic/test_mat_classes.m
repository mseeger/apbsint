% Test script for @Mat subclasses: They are compared to @MatDef for
% sparse matrices.
% Classes tested: @MatEye, @MatDiag, @MatSub, @MatContainer.

num_it = 100;
typStr = {'MatEye','MatDiag','MatSub','MatContainer'};

for it=1:num_it
  n = randi([50,200],1);
  B1arr = []; % For MatContainer
  B2comb = [];
  for typ = 1:4
    switch typ
     case 1
      B1 = ept.MatEye(n);
      smat = sparse(eye(n));
      B2 = ept.MatDef(smat);
      B2comb = [B2comb; smat];
      m=n;
     case 2
      v = 10*randn(n,1);
      B1 = ept.MatDiag(v);
      smat = sparse(diag(v));
      B2 = ept.MatDef(smat);
      B2comb = [B2comb; smat];
      m=n;
     case 3
      pm = randperm(n);
      m = randi([1,n],1);
      ind = pm(1:m)';
      B1 = ept.MatSub(n,ind);
      % Id_{ind,.}
      tmat = zeros(m,n);
      tmat(m*(ind-1)+(1:m)') = 1;
      smat = sparse(tmat);
      B2 = ept.MatDef(smat);
      B2comb = [B2comb; smat];
     case 4
      B1 = ept.MatContainer(B1arr);
      B2 = ept.MatDef(B2comb);
      m = size(B2comb,1);
    end
    if typ<4
      B1arr{typ} = B1;
    end
    % mvm
    v = randn(n,1);
    bv1 = mvm(B1,v);
    bv2 = mvm(B2,v);
    df = max(ept.reldiff(bv1,bv2));
    if df>1e-5
      fprintf(1,'mvm(%s): max(rdiff)=%f\n',typStr{typ},df);
    end
    v = randn(m,1);
    bv1 = mvm(B1',v);
    bv2 = mvm(B2',v);
    df = max(ept.reldiff(bv1,bv2));
    if df>1e-5
      fprintf(1,'mvmT(%s): max(rdiff)=%f\n',typStr{typ},df);
    end
    % getCol
    i = randi([1,n],1);
    bv1 = getCol(B1,i);
    bv2 = getCol(B2,i);
    df = max(ept.reldiff(bv1,bv2));
    if df>1e-5
      fprintf(1,'getCol(%s): max(rdiff)=%f\n',typStr{typ},df);
    end
    i = randi([1,m],1);
    bv1 = getCol(B1',i);
    bv2 = getCol(B2',i);
    df = max(ept.reldiff(bv1,bv2));
    if df>1e-5
      fprintf(1,'getColT(%s): max(rdiff)=%f\n',typStr{typ},df);
      %pause;
    end
    % matBtDgB
    v = randn(m,1);
    S1 = matBtDgB(B1,v);
    S2 = matBtDgB(B2,v);
    df = max(ept.reldiff(S1,S2));
    if df>1e-5
      fprintf(1,'matBtDgB(%s): max(rdiff)=%f\n',typStr{typ},df);
      %pause;
    end
    % diagBSymBt
    S = randn(n,n);
    bv1 = diagBSymBt(B1,S);
    bv2 = diagBSymBt(B2,S);
    df = max(ept.reldiff(bv1,bv2));
    if df>1e-5
      fprintf(1,'diagBSymBt(%s): max(rdiff)=%f\n',typStr{typ},df);
    end
  end
end
