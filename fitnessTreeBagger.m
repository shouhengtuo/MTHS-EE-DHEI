function accuracy = fitnessTreeBagger(X,Y)  
  model = TreeBagger(10,X,Y,'OOBPrediction','On', 'Method','classification');
  err = oobError(model);
  accuracy=1-err(end);
end



% X = data(:,88:100);
% Y = data(:,end);
% model = TreeBagger(30,X,Y,'OOBPrediction','On', 'Method','classification');
% % view(model.Trees{1},'Mode','graph')
% err = oobError(model);
% err(end)
% 
% 
% x = load('GeneExpression3.mat');
% data = x.diffData';
% label = x.labels;
% X = data(:,[1 35 202 129]);
% Y = label;
% model = TreeBagger(30,X,Y,'OOBPrediction','On', 'Method','classification');
% err = oobError(model);
% err(end)

