function [output] = nearneigh (MatA, dim)

sizA = size(MatA); output = zeros(sizA);
for k = 1:sizA(1)
    for l = 1:sizA(2)
        for m = 1:sizA(3)
            for o = 1:sizA(4)
                        
            on = MatA(k,l,m,o); 
            
            if dim == 1
                if k == sizA(1)
                    output(k,l,m,o) = abs(MatA(k-1,l,m,o)-on);
                elseif k == 1
                    output(k,l,m,o) = abs(MatA(k+1,l,m,o)-on);
                else
                up = MatA(k+1,l,m,o)-on; down = MatA(k-1,l,m,o)-on;
                 if abs(up)>=abs(down)
                output(k,l,m,o) = abs(down);
                 else 
                output(k,l,m,o) = abs(up);
                 end
                end
            elseif dim == 2
                 if l == sizA(2)
                    output(k,l,m,o) = abs(MatA(k,l-1,m,o)-on);
                 elseif l == 1
                    output(k,l,m,o) = abs(MatA(k,l+1,m,o)-on);
                 else
                up = MatA(k,l+1,m,o)-on; down = MatA(k,l-1,m,o)-on;
                  if abs(up)>=abs(down)
                output(k,l,m,o) = abs(down);
                 else 
                output(k,l,m,o) = abs(up);
                  end
                 end
            elseif dim == 3
                 if m == sizA(3)
                    output(k,l,m,o) = abs(MatA(k,l,m-1,o)-on);
                 elseif m == 1
                    output(k,l,m,o) = abs(MatA(k,l,m+1,o)-on);
                 else 
                up = MatA(k,l,m+1,o)-on; down = MatA(k,l,m-1,o)-on;
                 if abs(up)>=abs(down)
                output(k,l,m,o) = abs(down);
                 else 
                output(k,l,m,o) = abs(up);
                 end
                 end
            elseif dim == 4
                 if o == sizA(4)
                    output(k,l,m,o) = abs(MatA(k,l,m,o-1)-on);
                 elseif o == 1
                    output(k,l,m,o) = abs(MatA(k,l,m,o+1)-on);
                 else 
                up = MatA(k,l,m,o+1)-on; down = MatA(k,l,m,o-1)-on;
                 if abs(up)>=abs(down)
                output(k,l,m,o) = abs(down);
                 else 
                output(k,l,m,o) = abs(up);
                 end
                 end
                 
            end
                            
            end   
        end
    end
end
end

