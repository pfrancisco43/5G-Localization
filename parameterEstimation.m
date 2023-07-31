% Author: Paulo Francisco
% Objective: Estimation of the parameters
% Syntax:
%       [ToAs_E, AoDs_E, AoAs_E]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, H, fin)
% Inputs:
%       Nt, Nr, Ns: Antennas Tx, Rx and number of symbols
%       N, B, c: Numbert of subcariers, BW and propagation speed
%       y, x, L: signal transmited, original symbols and number of paths
%       L_az, L_el, H, fin: Length of Grid sensing, Channel model and option of fine tunning
%
% Outputs:
%       ToAs_E, AoDs_E, AoAs_E - Parameters estimated
%
function [ToAs_E, AoDs_E, AoAs_E]=parameterEstimation (Nt, Nr, Ns, N, B, c, y, x, L, L_az, L_el, H, fin)
    Ts=1/B;         % Sampling period in us

    array_Tx=getArray(Nt);
    array_Rx=getArray(Nr);
%Grid for the general scenario:    
%     rangeAz_aod=linspace(1,180,L_az)*pi/180;
%     rangeEl_aod=linspace(91,180,L_el)*pi/180;
%     rangeAz_aoa=[linspace(-89,0,L_az/2) linspace(181,270,L_az/2)]*pi/180;
%     rangeEl_aoa=linspace(1,90,L_el)*pi/180;

%Grid - In the specific simulation scenario (limit the grid for increase speed):
    rangeAz_aod=linspace(10,120,L_az)*pi/180;
    rangeEl_aod=linspace(91,120,L_el)*pi/180;
    rangeAz_aoa=[linspace(-60,-30,L_az/2) linspace(181,215,L_az/2)]*pi/180;
    rangeEl_aoa=linspace(58,85,L_el)*pi/180;

    %% Create dictionary 
    k=1;
    NrA=length(rangeAz_aoa);
    NrE=length(rangeEl_aoa);
    Ur=zeros(Nr,NrA*NrE);
    for i=1:NrA
        azi=rangeAz_aoa(i);
        for j=1:NrE
            ele=rangeEl_aoa(j);
            Ur(:,k)=getResponse(array_Rx,[azi ele])*sqrt(Nr);
            k=k+1;
        end
    end
    %

    k=1;
    Ut=zeros(Nt,NrA*NrE);
    for i=1:NrA
        azi=rangeAz_aod(i);
        for j=1:NrE
            ele=rangeEl_aod(j);
            Ut(:,k)=getResponse(array_Tx,[azi ele])*sqrt(Nt);
            k=k+1;
        end
    end


    %% Vectorize and generation of the basis
    yb=zeros(Nr*Ns,N);
    Omega=zeros(Nr*Ns,(NrA*NrE)^2,N);
    for n=1:N
        yb(:,n)=reshape(y(:,:,n),Nr*Ns,1); 
        Omega(:,:,n)=kronMult((Ut'*x(:,:,n)).',Ur);
    end


    %% Visualize the signal channel in AOA/AOD space
    HH=mean(H,3);
    spec=abs(Ur'*HH*Ut);
    mesh(abs(spec));
    xlabel('AOD'); ylabel('AOA');


    %% run DCS-SOMP
    [indices,h_hat]=DCSSOMP(yb,Omega,L); 
    
    %% compute the distance
    distances=zeros(1,L);
    for l=1:L
        distances(l)=-mean(diff(phase(h_hat(l,:))))*(N*Ts)*c/(2*pi); 
        if (distances(l)<0)
            distances(l)=distances(l)+N*Ts*c;
        end
    end

    %% estimate angles
    aodEst=zeros(L,2);
    aoaEst=zeros(L,2);
    for l=1:L
        id_1=ceil(indices(l)/(NrA*NrE)); 
        if mod(id_1,NrE)== 0
            id_A=id_1/NrE; 
            id_E=mod(id_1,NrE)+NrE;
        else
            id_A=fix(id_1/NrE+1);
            id_E=mod(id_1,NrE);
        end
        aod_az=rangeAz_aod(id_A); 
        aod_el=rangeEl_aod(id_E); 
        aodEst(l,1:2)=[aod_az aod_el];

        id_1=indices(l)-(id_1-1)*(NrA*NrE); 
        if mod(id_1,NrE)== 0
            id_A=id_1/NrE; 
            id_E=mod(id_1,NrE)+NrE; )
        else
            id_A=fix(id_1/NrE+1);
            id_E=mod(id_1,NrE);
        end
        aoa_az=rangeAz_aoa(id_A); 
        aoa_el=rangeEl_aoa(id_E); 
        aoaEst(l,1:2)=[aoa_az aoa_el];
    end

    %this part is optional and does not include in the paper
    %use the DCS-SOMP to refine results
    if fin==1
        for i=1:L
            [aodEst(i,:),aoaEst(i,:),h]=paramEstimateFining(aodEst(i,:),aoaEst(i,:),array_Tx,array_Rx,x,N,yb,Nt,Nr);
            distances(i)=-mean(diff(phase(h)))*(N*Ts)*c/(2*pi); 
            if (distances(i)<0)
                distances(i)=distances(i)+N*Ts*c;
            end
        end
    end

    %Order the parameters
    [ToAs_E, AoDs_E, AoAs_E]=orderParameters(distances'/c, aodEst, aoaEst);

end

function [ToAs_E, AoDs_E, AoAs_E]=orderParameters(toas, aods, aoas)
    [ToAs_E,i] = sort(toas);
    AoDs_E = aods(i,:);
    AoAs_E = aoas(i,:);
end


function X = kronMult(A,B)

    %KRON Kronecker product.
    %   kron(A,B) returns the Kronecker product of two matrices A and B, of
    %   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
    %   block matrix in which the (i,j)-th block is defined as A(i,j)*B.
    %   Version: 06/02/2011
    %   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
    [I, J] = size(A);
    [K, L] = size(B);
    if ~issparse(A) && ~issparse(B)

        % Both matrices are dense.
        A = reshape(A,[1 I 1 J]);
        B = reshape(B,[K 1 L 1]);
        X = reshape(bsxfun(@times,A,B),[I*K J*L]);

    else

        % One of the matrices is sparse.
        [ia,ja,sa] = find(A);
        [ib,jb,sb] = find(B);
        ix = bsxfun(@plus,K*(ia(:)-1).',ib(:));
        jx = bsxfun(@plus,L*(ja(:)-1).',jb(:));

        % The @and operator is slightly faster for logicals.
        if islogical(sa) && islogical(sb)
            X = sparse(ix,jx,bsxfun(@and,sb(:),sa(:).'),I*K,J*L);
        else
            X = sparse(ix,jx,double(sb(:))*double(sa(:).'),I*K,J*L);
        end
    end
end