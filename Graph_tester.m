
for ind =1:5:size(M100dat,2)
offset = 5
plot(M100dat(:, ind:ind+offset), 'Linewidth',2); legend(num2str((ind:ind+offset)'))
pause
end
%%
vector = [50, 91]
plot(M100dat(:,vector), 'Linewidth',2); legend(num2str((vector)'))