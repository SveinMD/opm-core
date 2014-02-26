function array = readXmlDataArray(fileName)
xml = xmlread(fileName);

dataArrays = xml.getElementsByTagName('DataArray');

for i = 0:dataArrays.getLength-1
    da = dataArrays.item(i);
    if(nameIs('saturation', da))
        array = readNode(da);
    end
end

end

function array = readNode(dataArray)
array = [];
contents = dataArray.getTextContent;
cells = textscan(char(contents),'%f %f');
array = [cells{1} cells{2}];
end

function nameMatch = nameIs(value, dataArray)
nameMatch = false;

if(dataArray.hasAttributes)
    attributes = dataArray.getAttributes;
    natt = attributes.getLength;
    for i = 0:natt-1
        att = attributes.item(i);
        if(strcmp(att.getName, 'Name') && strcmp(att.getValue, value))
            nameMatch = true;
            return;
        elseif(strcmp(att.getName, 'Name'))
            nameMatch = false;
            return;
        end
    end
end

end