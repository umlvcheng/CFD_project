function printbars(M)
  fprintf(1,'%s', '       |');
  for i = 2:M+1
    fprintf(1,'%s', '                   |');
  end
    fprintf(1,'\n');
