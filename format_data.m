function [Y_train,X_train,train_inds,train_time,nstates_X,nstates_Y,X_in,sigma_X,sigma_Y,nX,nY] = ...
    format_data(Ycell,X,t_epi,subsamp_t,train_time_frac,num_train_inds,test_length,snr_X,snr_Y,noise_alg_X,noise_alg_Y,seed1,seed2)

    if ~exist('seed1','var')
        seed1='shuffle';
    end
    if ~exist('seed2','var')
        seed2='shuffle';
    end
    if isempty(seed1)
        seed1=seed2;
    elseif length(seed1)>1
        num_train_inds = length(seed1);
    end

    nstates_X = size(X,2);
    nstates_Y = size(Ycell{1},2);
    test_length = max(min(test_length, length(Ycell)-1-num_train_inds),0);
    
    if or(length(seed1)==1,isequal(seed1,'shuffle'))
        rng(seed1)
        X_in = randperm(length(Ycell)-2-test_length,min(length(Ycell)-2-test_length,num_train_inds-1))';
        X_in = sort([1;X_in+1],'ascend');
    else
        X_in = seed1(:);
    end
    Y_train = arrayfun(@(n)Ycell{n}(1:subsamp_t:floor(end*train_time_frac),:),X_in,'uni',0);
    train_inds = unique([X_in(:);X_in(:)+1]);
    X_train = X(train_inds,:);
    % train_time = t_epi(1:end*train_time_frac);
    train_time = cellfun(@(t)t(1:subsamp_t:floor(end*train_time_frac)),t_epi(X_in),'uni',0);
    X_in = find(ismember(train_inds,X_in));

    rng(seed2)
    if isequal(noise_alg_X,'AWGN')
        sigma_X = snr_X*std(X_train);
        X_train = X_train+sigma_X.*randn(size(X_train));
    elseif isequal(noise_alg_X,'MWGN')
        sigma_X = snr_X + rms(X_train)*0;
        X_train = X_train.*(1+sigma_X*randn(size(X_train)));
    elseif isequal(noise_alg_X,'logn')
        % expsig = fzero(@(x)x.^4-2*x+1-snr_X^2,[1,2^(1/3)]);
        % sig = sqrt(log(expsig)*2);
        % sigma_X = sig+rms(X_train)*0;
        % X_train = X_train.*exp(sigma_X.*randn(size(X_train)));
        sigma_X = repmat(snr_X*rms(X_train),size(X_train,1),1);
        if snr_X>0
            X_train = X_train.^2./sqrt(X_train.^2+rms(X_train).^2*snr_X^2).*exp(sqrt(log(1+rms(X_train).^2./X_train.^2*snr_X^2)).*randn(size(X_train)));
        end
    end
    if isequal(noise_alg_Y,'AWGN')
        sigma_Y = cellfun(@(y)snr_Y*std(cell2mat(Y_train)),Y_train,'uni',0);
        % sigma_Y = cellfun(@(y)snr_Y*std(y),Y_train,'uni',0);
        Y_train = cellfun(@(y,s) y+s.*randn(size(y)),Y_train,sigma_Y,'uni',0);
    elseif isequal(noise_alg_Y,'MWGN')
        % sigma_Y = cellfun(@(y)snr_Y*std(cell2mat(Y_train)),Y_train,'uni',0);
        sigma_Y = cellfun(@(y)snr_Y+std(y)*0,Y_train,'uni',0);
        Y_train = cellfun(@(y,s) y+s.*randn(size(y)),Y_train,sigma_Y,'uni',0);
    elseif isequal(noise_alg_Y,'logn')
        % expsig = fzero(@(x)x.^4-2*x+1-snr_Y^2,[1,2^(1/3)]);
        % sig = sqrt(log(expsig)*2);
        % sigma_Y = cellfun(@(y)sig+std(y)*0,Y_train,'uni',0);
        % Y_train = cellfun(@(y,s) y.*exp(s.*randn(size(y))),Y_train,sigma_Y,'uni',0);
        sigma_Y = cellfun(@(y)snr_Y*rms(y),Y_train,'uni',0);
        if snr_Y>0
            Y_train = cellfun(@(y)y.^2./sqrt(y.^2+rms(y).^2*snr_Y^2).*exp(sqrt(log(1+rms(y).^2./y.^2*snr_Y^2)).*randn(size(y))),Y_train,'uni',0);
        end
        % Y_train = cellfun(@(y)y.*exp(-0.5*log(1+snr_Y^2)+sqrt(log(1+y*snr_Y^2)).*randn(size(y))))),Y_train,'uni',0);
    end

    nX = mean(abs(X_train));
    X_train = X_train./nX;
    nY = mean(abs(cell2mat(Y_train)));
    Y_train = cellfun(@(x)x./nY,Y_train,'uni',0);
    sigma_X = sigma_X./nX;
    sigma_Y = cellfun(@(s)num2cell(s./nY),sigma_Y,'uni',0);

end