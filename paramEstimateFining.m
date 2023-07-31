% Author: Paulo Francisco
% Objective: Fine Tunning to parameter estimation
% This function is in test

function [aodR, aoaR, h]=paramEstimateFining(aod,aoa,array_Tx,array_Rx,F,K,R,Nt,Nr)
    res=deg2rad(0.001);

    %improve azimuth AoD
    ang=aod(1);
    ele=aod(2);
    Ur=getResponse(array_Rx,[aoa(1) aoa(2)])*sqrt(Nr);
    maxval=0;
    while 1
        Ut=getResponse(array_Tx,[ang ele])*sqrt(Nt);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestAz=ang;
        else
            break;
        end
        ang=ang-res;
        if ang <= 0 || ang >=pi
            break
        end
    end
    ang=aod(1)+res;
    while 1
        Ut=getResponse(array_Tx,[ang ele])*sqrt(Nt);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestAz=ang;
        else
            break;
        end
        ang=ang+res;
        if ang <= 0 || ang >=pi
            break
        end
    end
    %#################################################

    %improve elevation AoD
    azi=bestAz;
    ang=aod(2);
    Ur=getResponse(array_Rx,[aoa(1) aoa(2)])*sqrt(Nr);
    maxval=0;
    while 1
        Ut=getResponse(array_Tx,[azi ang])*sqrt(Nt);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestEl=ang;
        else
            break;
        end
        ang=ang-res;
        if ang <= 0 || ang >=pi/2
            break
        end
    end
    ang=aod(2)+res;
    while 1
        Ut=getResponse(array_Tx,[azi ang])*sqrt(Nt);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestEl=ang;
        else
            break;
        end
        ang=ang+res;
        if ang <= 0 || ang >=pi/2
            break
        end
    end
    %#################################################
    aodR=[bestAz bestEl];


    %improve azimuth AoA
    ang=aoa(1);
    ele=aoa(2);
    Ut=getResponse(array_Tx,[aodR(1) aodR(2)])*sqrt(Nt);
    maxval=0;
    while 1
        Ur=getResponse(array_Rx,[ang ele])*sqrt(Nr);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestAz=ang;
        else
            break;
        end
        ang=ang-res;
        if ang <= 0 || ang >=pi
            break
        end
    end
    ang=aoa(1)+res;
    while 1
        Ur=getResponse(array_Rx,[ang ele])*sqrt(Nr);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestAz=ang;
        else
            break;
        end
        ang=ang+res;
        if ang <= 0 || ang >=pi
            break
        end
    end
    %#################################################

    %improve elevation AoA
    azi=bestAz;
    ang=aoa(2);
    Ut=getResponse(array_Tx,[aodR(1) aodR(2)])*sqrt(Nt);
    maxval=0;
    while 1
        Ur=getResponse(array_Rx,[azi ang])*sqrt(Nr);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestEl=ang;
        else
            break;
        end
        ang=ang-res;
        if ang <= 0 || ang >=pi/2
            break
        end
    end
    ang=aoa(2)+res;
    while 1
        Ur=getResponse(array_Rx,[azi ang])*sqrt(Nr);

        for k=1:K
            om(:,k)=kron((Ut'*F(:,:,k)).',Ur);
        end

        cost=0;
        for k=1:K
            cost=cost+abs(om(:,k)'*R(:,k))/norm(om(:,k));
        end

        if cost > maxval
            maxval=cost;
            bestEl=ang;
        else
            break;
        end
        ang=ang+res;
        if ang <= 0 || ang >=pi/2
            break
        end
    end
    %#################################################
    aoaR=[bestAz bestEl];

    %improve ToA
    Ut=getResponse(array_Tx,[aodR(1) aodR(2)])*sqrt(Nt);
    Ur=getResponse(array_Rx,[aoaR(1) aoaR(2)])*sqrt(Nr);

    h=zeros(1,k);
    for k=1:K
        psi=kron((Ut'*F(:,:,k)).',Ur);
        beta=psi'*R(:,k)/(norm(psi))^2;
        [Q,Rqr]=qr(psi,0);
%         h(k)=(Rqr\Q')*psi*beta;
        h(k)=inv(Rqr)*Q'*psi*beta;
    end

end