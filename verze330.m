clc; clear all; close all;

%% UnIT Image processing
%  autoøi: Cindy Veselá, Veronika Pošustová,
%          Michaela Hláèiková, Jana Schwarzerová
%                                              12.-13. 4. 2019
%                                              [Matlab2018a / Matlab2015a]

% Zadání:
% Aplikace, která dokáže naèíst obrázek  TEM ve formátu TIFF.
% Analyzovat zda je èi není vstupní obraz viditelný elektronový svazek a
% proložit jeho okraje elipsou.
% Vstup ... obrázek z TEM mikroskopie
% Výstup ... parametry elipsy

addpath('D:\4. roèník\UNIT\IP-EKO-Cider-Gang-master\unit-2019-test\data\test')
fileFolder = 'D:\4. roèník\UNIT\IP-EKO-Cider-Gang-master\unit-2019-test\data\test';
dirOutput = dir(fullfile(fileFolder,'*.tiff'));
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);

obr = cell(1,numFrames);
vektor = cell(numFrames+1,7);
%vektor(1,:) = ['filename','ellipse_centre_x','ellipse_centre_y','ellipse_majoraxis','ellipse_minoraxis','ellipse_angle','elapsed_time'];
%vektor(1,:) = [num2str(filename),num2str(ellipse_centre_x),num2str(ellipse_centre_y),num2str(ellipse_majoraxis),num2str(ellipse_minoraxis),num2str(ellipse_angle),num2str(elapsed_time)];
for i = 1:numFrames
    I = imread(fileNames{i});
    obr{1,i} = I;
end

obr2 = obr;
pom = cell2struct(obr,'data',1522);
puvodni_snimky = pom;

for i = 1:numFrames
    oo = pom(i).data;
    oo = im2double(oo);
    h = oo(:);
    maximum = max(h);
    minimum = min(h);
    novy_obr = (oo - minimum)/(maximum-minimum);
    puvodni_snimky(i).data = novy_obr;
end

clear oo
clear pom
clear maximum
clear minimum
clear h
clear obr2
clear novy_obr

%normalizace novy = (obr-min)/(max-min) (stejne jako imshow(obr,[])
    pocitadlo = 0;

for i = 1:numFrames
    tic
    obr = puvodni_snimky(1,i).data;
    T = graythresh(obr);
    binaryImage = im2bw(obr,T);
    binaryImage = imfill(binaryImage,'holes');
    
    s1 = regionprops(binaryImage,'Area','Centroid');
    s = regionprops(binaryImage,'Area');
   % pom_m = cell2mat( struct2cell( s ) );
    %pom_m = max(pom_m);
   
    
        
        m = cell2mat( struct2cell( s ) );
        m = max(m);
        prah = round(0.6*m);
        BW = bwareaopen(binaryImage,prah);
        s = regionprops(BW,{...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation'});
        
        % ovìøení kulovitosti:
        [B2,L2] = bwboundaries (BW,'noholes');
        labeledImage8 = bwlabel(BW, 8);
        labeledImage = labeledImage8;
        BW2 = labeledImage >= 1;
        binaryImage = imclearborder(BW2);
        binaryImage = bwareafilt(BW2, 3);
        [labeledImage, numberOfBlobs] = bwlabel(binaryImage, 8);
        props = regionprops(binaryImage, 'Area', 'Perimeter');
        allAreas = [props.Area];
        allPerimeters = [props.Perimeter];
        circularities = allPerimeters .^ 2 ./ (4 * pi * allAreas);
        maxAllowableCircularity = 6; %kulatý objekt
        keeperIndexes = find(circularities <= maxAllowableCircularity);
        binaryImage = ismember(labeledImage, keeperIndexes);
        binaryImage(binaryImage>0)=1;
        binaryImage = imfill(binaryImage,'holes');
        %%spojení max objektu s nekruhovitìjším objektem
        s = regionprops(binaryImage,{...
            'Extrema',...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation'});
        pocet = size(s);
        
        if pocet(1) > 1 
            obr = puvodni_snimky(1,i).data;
            T = graythresh(obr);
            binaryImage = im2bw(obr,T);
            binaryImage = imfill(binaryImage,'holes');
            s1 = regionprops(binaryImage,'Area','Centroid');
            s = regionprops(binaryImage,'Area');
            pom_m = cell2mat( struct2cell( s ) );
            pom_m = min(pom_m);
            
           % vyuzij jiny algoritmus - janca
            if length(s)> 3 && pom_m > 250
                pom = [s1(:).Centroid];
                new_Centroid(1) = sum(pom(1:2:end))/(length([s1(:).Centroid])/2);
                new_Centroid(2) = sum(pom(2:2:end))/(length([s1(:).Centroid])/2);
                
                [rad_pom,sloup_pom] = find(binaryImage == 1);
                pom_rad_min = min(rad_pom);
                hled_pix_rad(1) = rad_pom(min(find(rad_pom==pom_rad_min)));
                hled_pix_rad(2) = sloup_pom(min(find(rad_pom==pom_rad_min)));
                
                pom_sloup_min = min(sloup_pom);
                hled_pix_sloup(1) = rad_pom(min(find(sloup_pom==pom_sloup_min)));
                hled_pix_sloup(2) = sloup_pom(min(find(sloup_pom==pom_sloup_min)));
                
                prumer = sqrt((hled_pix_sloup(1)-hled_pix_sloup(2))^2+(new_Centroid(1)-new_Centroid(2))^2);
                cas = toc*1000;
                vektor(i,:) = [fileNames(i),num2str(new_Centroid(1)),num2str(new_Centroid(2)),num2str(prumer/2),num2str(prumer/2),num2str(0),num2str(cas)];
                pocitadlo = pocitadlo + 1;
            else
                 rozmer = size(obr);
                 cas = toc*1000;
                 vektor(i,:) = [fileNames(i),' ',' ',' ',' ',' ',num2str(cas)];
            end 
           
%             disp('moc centroidu')
%             rozmer = size(obr);
%             vektor(i,:) = [fileNames(i),' ',' ',' ',' ',' ',num2str(rozmer(1)),num2str(rozmer(2))];
           
        elseif pocet(1) == 0
           % disp('zadny centroid')
            rozmer = size(obr);
            cas = toc*1000;
            vektor(i,:) = [fileNames(i),' ',' ',' ',' ',' ',num2str(cas)];
        else
         [radky,sloupce] = size(binaryImage);
         novy = zeros(radky,sloupce);
         novy = [binaryImage, novy];
         
         [z,a,b,alpha] = fitellipse([s.Extrema(:,1),s.Extrema(:,2)]);
         if alpha <0
         alpha = abs(radtodeg(alpha));
       
         end
         cas = toc*1000;
         vektor(i,:) = [fileNames(i),num2str(z(1)),num2str(z(2)),num2str(a/2),num2str(b/2),num2str(alpha), num2str(cas)];
         
        end
    end


tabulka = cell2table(vektor);
writetable(tabulka,'vysledky.csv');


