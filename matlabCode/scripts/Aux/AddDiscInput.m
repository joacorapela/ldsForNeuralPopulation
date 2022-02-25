function seq = AddDiscInput(seq)


for T=1:length(seq)
  Sum = (seq(T).u(1,:))+(seq(T).u(2,:));
  Dif = (seq(T).u(1,:))-(seq(T).u(2,:));
  
  seq(T).u(1:2,:) = [Sum;Dif];
    
end