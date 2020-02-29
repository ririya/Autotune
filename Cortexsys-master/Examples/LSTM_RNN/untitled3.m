indNumber = findstr(Xlong, 2);
xNumber = Xlong(indNumber');
xNumber= xNumber';

indBad1 = find(indNumber >= length(Xlong));
indBad2 = find(indNumber == 1);
indBad = [indBad1 indBad2];

indNumber(indBad) = [];

xNumber(indBad) = [];

indBad = find(indNumber == 1);

xPrevNumber = Xlong(indNumber-1);
xPrevNumber = xPrevNumber';


xNextNumber = Xlong(indNumber+1);
xNextNumber = xNextNumber';
xFinal = [xPrevNumber xNumber xNextNumber];

uniqueRows = unique(xFinal,'rows');

countRow = [];
for i=1:size(uniqueRows,1)
    countRow(i) = length(findstr(Xlong, uniqueRows(i,:)));
    
end

[uniqueRows countRow']
