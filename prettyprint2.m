function prettyprint2(fu,fp,c)

  M = length(fp);
  if (nargin > 2)
    for i = 1:M
      fprintf(1,' %s %7.2f %s %7.2f', c,fu(i),c,fp(i));
    end
      fprintf(1,' %s %7.2f %s\n', c, fu(M+1), c);
  else
    for i = 1:M
      fprintf(1,'%8.2f %8.2f', fu(i),fp(i));
    end
      fprintf(1,'%8.2f\n', fu(M+1));
  end
