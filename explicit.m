%   Copyright (C) 2020  Ma Shisong Lab  All Rights Reserved.

classdef explicit
	properties
		beta =  [];   % tf (row)  x gene (col)
		beta_pvalue = [];  % tf x gene
		TF_name = []; % the name of TF genes, plus intercept
		Target_name = []; % the name of target genes
		NRMSE = []; % normalized root mean squre error for the training dataset
		Correlation_by_sample = []; % between predicted gene expression values vs. actual valuses of every sample
		Correlation_by_target_gene = []; % between predicted values va. actual values of every target gene across all samples
		SST = [];   % by target gene
		SSR = [];   % by target gene
		SSE = [];   % by tareget gene
		Fstat = []; % by tareget gene
		Fpvalue = []; % by target gene
		SigEdges = []; % sinificant edges with pValue smaller than cutoff
	end

	methods (Static)
		function obj=explicit (x1, x2, tf_name, target_name) % x1, tf gene expression matrix; x2, target gene matrix. samples in rows, tfs/genes in column.
			B = x2;
			A = [ones(size(x1,1),1) x1];
			beta = (A' * A) \ (A' * B);

			Bp = A * beta;
			c = arrayfun(@(k) corr(Bp(k,:)',B(k,:)'),1:size(B,1),'Uni',1);
			obj.Correlation_by_sample = c';

			cc = arrayfun(@(k) corr(Bp(:,k),B(:,k)),1:size(B,2),'Uni',1);
			obj.Correlation_by_target_gene = cc';

			r = Bp - B;
			rt = r';
			Bt = B';
			obj.NRMSE = sqrt( sum(sum(rt.^2)) / sum(sum(Bt.^2)));
			clear Bp Bt r rt;

			H = A * inv( A' * A ) * A';
			I = eye(size(A, 1));
			J = ones(size(A,1)) / size(A,1);

			SST = diag( B' * ( I - J ) * B );
			SSR = diag( B' * ( H - J ) * B );
			SSE = diag( B' * ( I - H ) * B );

			clear H I J;

			dof_sse = size(A, 1) - size(A, 2);
			dof_ssr = size(A, 2) - 1;
			se_beta = sqrt( diag( inv(A' * A) ) * (SSE/dof_sse)');
			tStat = beta ./ se_beta;
			beta_pvalue = tcdf( -abs(tStat), dof_sse) * 2;
			Fstat = (SSR / dof_ssr) ./ (SSE/dof_sse);
			obj.Fstat = Fstat;
			obj.Fpvalue = fcdf(Fstat,dof_ssr,dof_sse,'upper');
			obj.beta = beta;
			obj.beta_pvalue = beta_pvalue;
			obj.SST = SST;
			obj.SSR = SSR;
			obj.SSE = SSE;

			idx = find(beta_pvalue <= 0.00001);
			e3 = beta(idx);
			e3 = round(e3,4);
			e4 = beta_pvalue(idx);
			e4 = round(e4,4,'significant');
			[e1, e2] = ind2sub(size(beta), idx);
			
			if nargin > 2
				tf_name = string(tf_name);
				tf_name = ["intercept" tf_name'];
				tf_name = regexp( tf_name,'([a-zA-Z0-9_\.]+)','once','match');
				target_name = string(target_name);
				target_name = target_name';
				target_name = regexp( target_name,'([a-zA-Z0-9_\.]+)','once','match');
				e1 = tf_name(e1);
				e2 = target_name(e2);
				e1 = e1';
				e2 = e2';
			end

			obj.TF_name = tf_name;
			obj.Target_name = target_name;
			colName = {'Gene';'TF';'beta';'beta_pvalue'};
			SigEdges = table(e2,e1,e3,e4,'VariableNames',colName);
			obj.SigEdges = SigEdges;

		end
	end
end

