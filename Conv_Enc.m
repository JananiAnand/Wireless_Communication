function Conv_Enc(image_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convolutional Encoding with Viterbi Decoding
%
%This program performs Convolutional Encoding with Viterbi Decoding
%It also evaluates the effect of noise of encoding
%
%function Conv_Enc(image_file)
%Eg: Conv_Enc('image1.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exist image_file;

if ~ans
    fprintf('\n')
    error('Run Conv_Enc(''image1.jpg'') or type help Conv_Enc for documentation')
end

%SNR Values
SNR = [0 3 6 9 12 15];
k = 0;
time = 5;
%Reading the image file
I_in = imread(image_file);
%Reshaping image matrix for processing
B = dec2bin(I_in);
C = reshape(B',1,numel(B));
data = C - '0';
%Trellis Generation
trellis = poly2trellis(4,{'1 + x + x^3 ','x + x^2','x^2 + x^3'});
%Convolutional Encoding
codedData = convenc(data,trellis);
%Viterbi Decoding
decodedData = vitdec(codedData,trellis,length(data),'trunc','hard');
%Number of Bit error between original and viterbi decoded image
x = biterr(data,decodedData);

fprintf('Q1. Convolutional Encoded Data \n')
fprintf('codedData = convenc(data,trellis)\n')
fprintf('\n')
pause(time)

fprintf('Q3. Viterbi Decoding\n')
fprintf('No of bit errors: %d\n',x)
fprintf('Original Image: Fig 1\n')
fprintf('Viterbi Decoded Image: Fig 2\n')
fprintf('\n')

decodedDataPlot = num2str(decodedData,'%d');
D = reshape(decodedDataPlot,size(B,2),size(B,1));
I_out = uint8(reshape(bin2dec(D'),size(I_in)));

%Displaying Original Image
figure('Name','Original Image','NumberTitle','off'); 
imshow(I_in);
title('Original Image')
pause(time)
%Displaying Decoded Image
figure('Name','Viterbi Decoded Image','NumberTitle','off');
imshow(I_out);
title('Viterbi Decoded Image')
pause(time)

codedDataBipolar = codedData;

%Biploar Signalling
for i = 1:length(codedDataBipolar)
    if codedDataBipolar(i) == 0
        codedDataBipolar(i) = -1;
    end
end


figure('Name','Decoded image for various SNR of convolutional encoded bits','NumberTitle','off'); 
fprintf('Q4.\n')
fprintf('Signal to Noise Ratio (SNR) is useful in judging the impact of noise on system performance\n')
fprintf('It is seen that image quality degrades with lower SNR\n')
fprintf('While hard decoding has a lower accuracy as compared to soft decoding\n')
fprintf('\n')
fprintf('Decoded image for various SNR of convolutional encoded bits:Fig 3\n')
fprintf('\n')
fprintf('Convolutional Encoded Bits\n')
%Adding White Gaussian Noise
for i = 0:3:15
    a=awgn(codedDataBipolar,i);
    %Hard Decision
    for j = 1:length(a)
        if a(j) < 0.9
            a(j) = 0;
        else
            a(j) = 1;
        end
    end
	
k = k+1;
%Viterbi Decoding
decodedDataBipolar = vitdec(a,trellis,length(data),'trunc','hard');
fprintf('No of Error Bits for %d dB = ',i);
 
biterror(k) = biterr(data,decodedDataBipolar);
fprintf('%d\n',biterror(k))

decodedDataBipolarPlot = num2str(decodedDataBipolar,'%d');
D = reshape(decodedDataBipolarPlot,size(B,2),size(B,1));
I_out = uint8(reshape(bin2dec(D'),size(I_in)));

subplot(6,1,k)
imshow(I_out)
title(['SNR ' num2str(i) ' dB'])
end
pause(time)
fprintf('\n')

dataBipolar = data;

%Biploar Signalling
for i = 1:length(dataBipolar)
    if dataBipolar(i) == 0
        dataBipolar(i) = -1;
    end
end

k = 0;


figure('Name','Decoded image for various SNR of original bits','NumberTitle','off'); 
fprintf('Q5.\n')
fprintf('Decoded image for various SNR of original bits:Fig 4\n')
fprintf('\n')
fprintf('Original Bits\n')

%Adding White Gaussian Noise
for i = 0:3:15
    a=awgn(dataBipolar,i);
    %Hard Decision
    for j = 1:length(a)
        if a(j) < 0.9
            a(j) = 0;
        else
            a(j) = 1;
        end
    end
k = k+1;


fprintf('No of Error Bits for %d dB = ',i);
biterrorOrig(k) = biterr(data,a);
fprintf('%d\n',biterrorOrig(k))

dataBipolarPlot = num2str(a,'%d');
D = reshape(dataBipolarPlot,size(B,2),size(B,1));
I_out = uint8(reshape(bin2dec(D'),size(I_in)));

subplot(6,1,k)
imshow(I_out)
title(['SNR ' num2str(i) ' dB'])
end
pause(time)

fprintf('\n')
fprintf('BER vs SNR for convolutional encoded bits and original bits: Fig 5\n')
fprintf('\n')
fprintf('It is seen that image quality degrades severely with lower SNR when the original bits are taken as compared to the encoded bits\n')
figure('Name','BER vs SNR for convolutional encoded bitsand original bits','NumberTitle','off');
set(gca,'XLim',[1 10])
set(gca,'XTick',(1:1:10))
semilogy(biterror,10.^SNR,'-*')
xlabel('No of Bit Errors')
ylabel('SNR in dB')
hold on
semilogy(biterrorOrig,10.^SNR,'-*')
legend('No of Bit Error using Vertibri Decoding','No of Bit Error using Original Data')
hold off       