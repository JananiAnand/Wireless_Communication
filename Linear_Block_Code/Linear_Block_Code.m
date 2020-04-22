%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linear_Block_Code.m
%
%This script is used to implement Linear Block Codes with interleaving and framing
%This code uses coding parameters (n,k,d) to perform linear block coding.
%It is used to correct upto single bit errors
%
%Author: Janani
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Linear Block Code\n')
fprintf('\n')
fprintf('Coding Parameters:\n')
fprintf('\tn -> Code word length\n')
fprintf('\tk ->  Number of message bits\n')
fprintf('\td -> Hamming Distance\n')
fprintf('\tWe have used (n,k,d) -? (15,11,3)\n')
fprintf('\tA simple linear block code implementation can detect and correct\n')
fprintf('\tonly single error in the bit stream. If we have multi – bit errors \n')
fprintf('\tin burst, then it would lead to erroneous coding of bits.\n')
fprintf('\n')    

%Pause Time
time = 1;
pause(time);

%n bit codeword
n = 15;
%k message bits
k = 11;
%Hamming distance
hd = 3;
%k x k Identity Matrix
Ik = eye(k);
%(n-k) x (n-k) Identity Matrix
Ink = eye(n-k);
%n x n Identity Matrix for error pattern
E = eye(n);
%k x (n-k) A matrix
%Created from parity equation
%P1 = D1 + D4 + D5 + D7 + D9 + D10 + D11
%P2 = D1 + D2 + D4 + D6 + D7 + D8 + D9
%P3 = D2 + D3 + D5 + D7 + D8 + D9 + D10
%P4 = D3 + D4 + D6 + D8 + D9 + D10 + D11
A = transpose([1 0 0 1 1 0 1 0 1 1 1;1 1 0 1 0 1 1 1 1 0 0;0 1 1 0 1 0 1 1 1 1 0;0 0 1 1 0 1 0 1 1 1 1]);

fprintf('Generator Matrix\n')
fprintf('\tG = [I | A]\n')
fprintf('\tG is a k x n matrix\n')
fprintf('\tI is an k x k identity matrix\n')
fprintf('\tA is a k x (n-k) matrix\n')

pause(time);

%k x n Generator Matrix
%G = [Ik | A]
G = [Ik A];

fprintf('G =\n')
disp(G)
fprintf('\n')
fprintf('A =\n')
disp(A)
pause(time)
fprintf('\n')
fprintf('Channel Encoding\n')
fprintf('\n')
fprintf('\tParity Check Matrix\n')

%(n-k) x n Parity Check Matrix
H = [A.' Ink];

fprintf('H =\n')
disp(H)
fprintf('\n')

%Message Input
msg = 'EECS-202 - Basic Digital Communication with Networking';
%No. of binary bits representation in each ASCII value
nbits = 7;
%Message bits
D = [];
%Cell array for Syndrome matrix
S = {};
%Cell array for Codeword generation
C = {};
%Syndrome of received word
S_recv = {};
%Received message
R = {};
%Decoded ASCII equivalent
decode = [];
j = 0;
X = [];
R2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                CHANNEL ENCODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\t Syndromes')
fprintf('\n')

%Generation of Syndrome Matrix
for i = 1:n
    S{i} = H * transpose(E(i,:));
end

fprintf('S(1)=\n')
disp(S{1}) 
fprintf('S(2)=\n')
disp(S{2})
fprintf('S(3)=\n')
disp(S{3})
fprintf('S(4)=\n')
disp(S{4})
fprintf('S(5)=\n')
disp(S{5})
fprintf('S(6)=\n')
disp(S{6})
fprintf('S(7)=\n')
disp(S{7})
fprintf('S(8)=\n')
disp(S{8})
fprintf('S(9)=\n')
disp(S{9})
fprintf('S(10)=\n')
disp(S{10})
fprintf('S(11)=\n')
disp(S{11})
fprintf('S(12)=\n')
disp(S{12})
fprintf('S(13)=\n')
disp(S{13})
fprintf('S(14)=\n')
disp(S{14})
fprintf('S(15)=\n')
disp(S{15})
fprintf('\n')

%Converting ASCII value of message bits into equivalent binary
for i = 1:length(msg)
   D = [D (dec2bin(double(msg(i)),nbits) - '0')];
end

%Length of message in binary
len_msg_bin = length(D);
%Index for excess message bits 
index_excess = k*(round(len_msg_bin/k));
%Excess Message
msg_excess = [D(index_excess+1:end) zeros(1,k-(len_msg_bin-index_excess))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                CODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Codeword generation
for i = 1:k:index_excess
    j = j + 1;
    C{j} = mod(D(i:i+10)*G,2); 
end

%Codeword generation for excess message bits
C{j+1} = mod(msg_excess(1,:)*G,2);

%Syndrome decoding 
for i = 1:35
	%Received message with Errors
    R{i} = mod(C{i}+E(1,:),2);
    S_recv{i} = mod(H*R{i}.',2);
	
    for j = 1:n
        a(j) = isequal(S_recv{i},S{j});
    end
    
    b = find(a);
	%Correcting the errors
    R1 = mod((R{i} + E(b,:)),2);
    X = horzcat(X,R1(1:11));
    
end

%Conversion into ASCII equivalent
for k = 1:7:385
    decode = horzcat(decode,char(bin2dec(num2str(X(k:k+7-1)))));
end

pause(time);

fprintf('Coding\n')
fprintf('\n')
fprintf('\tWe have 54 characters in the input text. We have represented each\n')
fprintf('\tcharacter with 7 bits. Hence we are expecting to get 54 * 7 = 378 bits.\n')
fprintf('\n')
%Dispaying the decoded output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code performs Linear Block Coding with Interleaving and Framing.
%It is used to correct multi - bit errors which occurs in burst.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%n bit codeword
n = 15;
%k message bits
k = 11;
%Hamming distance
hd = 3;
%k x k Identity Matrix
Ik = eye(k);
%(n-k) x (n-k) Identity Matrix
Ink = eye(n-k);
%n x n Identity Matrix for error pattern
E = eye(n);
%k x (n-k) A matrix
%Created from parity equation
%P1 = D1 + D4 + D5 + D7 + D9 + D10 + D11
%P2 = D1 + D2 + D4 + D6 + D7 + D8 + D9
%P3 = D2 + D3 + D5 + D7 + D8 + D9 + D10
%P4 = D3 + D4 + D6 + D8 + D9 + D10 + D11
A = transpose([1 0 0 1 1 0 1 0 1 1 1;1 1 0 1 0 1 1 1 1 0 0;0 1 1 0 1 0 1 1 1 1 0;0 0 1 1 0 1 0 1 1 1 1]);
%k x n Generator Matrix
%G = [Ik | A]
G = [Ik A];
%(n-k) x n Parity Check Matrix
H = [A.' Ink];
%Message Input
msg = 'EECS-202 - Basic Digital Communication with Networking';
%Minimum range for random number generation
random_no_min = 0;
%Maximum range for random number generation
random_no_max = 100;
%Minimum range for random bit generation
random_bits_min = -1;
%Maximum range for random bit generation
random_bits_max = 1;
%Random number generation
random_no = ceil((random_no_max-random_no_min).*rand(1,1) + random_no_min)
%Random bit generation
random_bits = ceil((random_bits_max-random_bits_min).*rand(1,random_no) + random_bits_min);
%Synchronization pattern
sync = '01111110' - '0';
%To identify frame start
N = length(sync) + random_no;
%No. of binary bits representation in each ASCII value
nbits = 7;
%Message bits
D = [];
%Cell array for Syndrome matrix
S = {};
%Cell array for Codeword generation
C = {};
%Cell array for Received word
Intrlv_wrd = {};
%No of Interleaving blocks
block = 0;
%Syndrome of received word
S_recv = {};
%Decoded ASCII equivalent
decode = [];
j = 0;
X = [];
I = [];
R = {};
Framed_bits = [];
y = [];
%Converting ASCII value of message bits into equivalent binary
for i = 1:length(msg)
   D = [D (dec2bin(double(msg(i)),nbits) - '0')];
end

%Length of message in binary
len_msg_bin = length(D);
%Index for excess message bits 
index_excess = k*(round(len_msg_bin/k));
%Excess Message
msg_excess = [D(index_excess+1:end) zeros(1,k-(len_msg_bin-index_excess))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                CHANNEL ENCODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generation of Syndrome Matrix
for i = 1:n
    S{i} = H * transpose(E(i,:));
end

%Codeword generation
for i = 1:k:index_excess
    j = j + 1;
    C{j} = mod(D(i:i+10)*G,2); 
end

%Codeword generation for excess message bits
C{j+1} = mod(msg_excess(1,:)*G,2);

pause(time)
fprintf('Interleaving\n')
fprintf('\n')
fprintf('\tInterleaving with a block size of 16 x 15 = 240 bits\n')
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                INTERLEAVING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interleaving with block size 16 x 15 = 240 bits
for i = 1:16:32
    %No. of blocks
    block = block + 1;
    %16 x 15 block 
    Ib = vertcat(C{i:i+15});
    disp(Ib)
    %Interleaved block
    I = Ib(1:end);
	
    %Applying Synchronisation bits along with random bits
    R{block} = [sync random_bits I];
    R_x = R{block};
    Framed_bits = horzcat(Framed_bits,R_x);
    Intrlv_wrd{block} = [R_x(1:234) ~R_x(235:249) R_x(250:end)];
	R_y = Intrlv_wrd{block};
	y = horzcat(y,R_y);
	
end

%Interleaving for excess message bits
Ib = vertcat(C{33:35});
Ib_excess = vertcat(Ib,zeros(13,15));
I = Ib_excess(1:end);
disp(Ib_excess)
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FRAMING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_excess = [sync random_bits I];

pause(time)
fprintf('Framing\n')
fprintf('\n')
Framed_bits = horzcat(Framed_bits,R_excess)
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                INTRODUCING ERRORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Intrlv_wrd{block+1} = [R_excess(1:234) ~R_excess(235:249) R_excess(250:end)];
R_y = Intrlv_wrd{block + 1};

pause(time)
fprintf('Introduce Errors\n')
fprintf('\n')
y = horzcat(y,R_y)
fprintf('\n')

pause(time)
fprintf('De-Framing\n')
fprintf('\tfunction de_framed_word = de_framing(Framed_word)\n')
fprintf('\tde_framed_word = Framed_word((length(sync) + random_no):end)\n')
fprintf('\tend\n')
fprintf('\n')

pause(time)
fprintf('De-Interleaving\n')
fprintf('\n')
fprintf('\tWe could correct all the consecutive bit errors. Bits from a particular\n') 
fprintf('\tcode word are transmitted sequentially, so a B bit burst produces multi-bit\n')
fprintf('\terrors (B is 15 in our case). A solution to this problem is to interleave\n') 
fprintf('\tbits from B different code words. Now a B- bit burst produces 1-bit errors\n') 
fprintf('\tin B different code words.\n')
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 					DE-FRAMING | DE-INTERLEAVING | DE-CODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Frame start detection and de-interleaving
for l = 1:3
    Recvd_dec = Intrlv_wrd{l};
    Recvd_dec1 = Recvd_dec(N+1:end);
    Recvd_wrd = reshape(Recvd_dec1,[16,15]);

    %Syndrome Decoding
    for i=1:16
        S_recv{i} = mod(H*Recvd_wrd(i,:).',2);
		
        for j = 1:n
            a(j) = isequal(S_recv{i},S{j});
        end
        
        b = find(a);
        
        if isempty(E(b,:))
            R2 = Recvd_wrd(i,:);
        else
			%Correcting the Errors
            R2 = mod((Recvd_wrd(i,:) + E(b,:)),2);
        end
        X = horzcat(X,R2(1:11));
    end
end
%Conversion into ASCII equivalent
for k = 1:nbits:525
    decode = horzcat(decode,char(bin2dec(num2str(X(k:k+7-1)))));
end

pause(time)
fprintf('De - coding\n')
fprintf('\n')
fprintf('\tYes, we were able to recover the entire message. Although we\n') 
fprintf('\thad to adjust the k bit frames of message blocks. On applying interleaving\n')
fprintf('\ttechnique we got only a single bit error in each code word. So, we could\n') 
fprintf('\tcorrect it using syndrome decoding\n')
fprintf('\n')
%Displaying the decoded output
fprintf('\t\t\t\t\t\t\t\tDecoded output\n\t\t\t')
disp(decode)
fprintf('\n')
