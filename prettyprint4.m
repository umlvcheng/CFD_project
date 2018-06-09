function prettyprint4(fz,fv,fu,fp,string)
  disp(string)
  [M,N] = size(fp);
  printbars(M);
  prettyprint2(fz(:,N+1),fv(:,N+1),'-');
  for j = N:-1:1
    prettyprint2(fu(:,j),fp(:,j),' ');
    prettyprint2(fz(:,j),fv(:,j),'-');
  end
  printbars(M);
