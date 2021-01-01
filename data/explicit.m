%   Copyright (C) 2020  Ma Shisong Lab  All Rights Reserved.

classdef explicit
	properties
		beta =  [];   % tf (row)  x gene (col)
		beta_pvalue = [];  % tf x gene
		NRMSE = [];    % RMSE  RMST  NRMSE
		NRMSE_all = [];
		Correlation_by_sample = []; % between predicted values vs. actual valuses of every experiment
		Correlation_by_gene = [];
		SST = [];   % by gene
		SSR = [];   % by gene
		SSE = [];   % by gene
		Fstat = []; % by gene
		Fpvalue = []; % by gene
		SigEdges = []; % edges with pValue smaller than cutoff;
	end

	methods (Static)
		function obj=explicit (x1, x2, tf_name, gene_name) % x1, tf matrix, x2, target gene matrix. sample in rows, tf/gene in column. tf_name and gene_name are optional.
			B = x2;
			A = [ones(size(x1,1),1) x1];
			beta = (A' * A) \ (A' * B);

			Bp = A * beta;
			c = arrayfun(@(k) corr(Bp(k,:)',B(k,:)'),1:size(B,1),'Uni',1);
			obj.Correlation_by_sample = c';

			cc = arrayfun(@(k) corr(Bp(:,k),B(:,k)),1:size(B,2),'Uni',1);
			obj.Correlation_by_gene = cc';

			r = Bp - B;
			Bp = [];
			rt = r';
			r = [];
			Bt = B';
			fn1 = sqrt(sum(rt.^2));
			fn2 = sqrt(sum(Bt.^2));
			fn3 = fn1 ./ fn2;
			n = [fn1; fn2; fn3];
			obj.NRMSE = n';
			obj.NRMSE_all = sqrt( sum(sum(rt.^2)) / sum(sum(Bt.^2)));
			Bt = [];
			rt = [];

			H = A * inv( A' * A ) * A';
			I = eye(size(A, 1));
			J = ones(size(A,1)) / size(A,1);

			SST = diag( B' * ( I - J ) * B );
			SSR = diag( B' * ( H - J ) * B );
			SSE = diag( B' * ( I - H ) * B );

			H = [];
			I = [];
			J = [];

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

			cAB = corr(A, B);
			e5 = cAB(sub2ind(size(cAB),e1,e2));
			e5 = round(e5, 4);
			
			if nargin > 2
				tf_name = string(tf_name);
				tf_name = ["intercept" tf_name'];
				tf_name = regexp( tf_name,'([a-zA-Z0-9_\.]+)','once','match');
				gene_name = string(gene_name);
				gene_name = gene_name';
				gene_name = regexp( gene_name,'([a-zA-Z0-9_\.]+)','once','match');
				e1 = tf_name(e1);
				e2 = gene_name(e2);
				e1 = e1';
				e2 = e2';
			end

			colName = {'Gene';'TF';'beta';'beta_pvalue';'correlation'};
			SigEdges = table(e2,e1,e3,e4,e5,'VariableNames',colName);
			obj.SigEdges = SigEdges;

		end
	end
end

			
