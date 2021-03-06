classdef explicit_kfcv
	properties
		Total_repeats = [];
		Fold_number = [];
		Correlation_of_target_gene = []; %between predicted expression values and actual values of every target genes (in rows) across test samples, in all CV runs (in col).
		Target_gene_name = [];
		AllEdges = []; % All edges (TF-target gene pairs) of the preidictor model
		CV_Stat = []; % Statistics of the CV runs
	end

	methods (Static)
		function obj = explicit_kfcv (x1, x2, tf_name, gene_name, x3, x4) % x1: tf gene expression matrix; x2: target gene expression matrix; x3: number of repeats; x4: number of folds
	
			fprintf('Building predictor model with full dataset...\n');	

			cross_repeat_num = 5;
			fold = 10;

			obj.Target_gene_name = gene_name;
			obj.Total_repeats = cross_repeat_num;
			obj.Fold_number = fold;

			if nargin == 6
				cross_repeat_num = x3;
				fold = x4;
			end

			obj.Total_repeats = cross_repeat_num;
			obj.Fold_number = fold;
			
			B = x2;
			A = [ones(size(x1,1),1) x1];
			beta = (A' * A) \ (A' * B);

			H = A * inv( A' * A ) * A';
			I = eye(size(A, 1));
			J = ones(size(A,1)) / size(A,1);

			SST = diag( B' * ( I - J ) * B );
			SSR = diag( B' * ( H - J ) * B );
			SSE = diag( B' * ( I - H ) * B );

			clear H I J;

			dof_sse = size(A, 1) - size(A, 2);
			se_beta = sqrt( diag( inv(A' * A) ) * (SSE/dof_sse)');
			tStat = beta ./ se_beta;
			beta_pvalue = tcdf( -abs(tStat), dof_sse) * 2;

			idx = find(beta_pvalue <= 2);
			e3 = beta(idx);
			e3 = round(e3,5);
			e4 = beta_pvalue(idx);
			e4 = round(e4,4,'significant');
			[e1, e2] = ind2sub(size(beta), idx);
				
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

			fprintf('\nStarting %d repeats of %d-fold cross-validation runs ...\n\n Rep-Iteration#\tNRMSE_Training\tNRMSE_Test\tR_Training\tR_Test', cross_repeat_num, fold);
			crv = zeros( cross_repeat_num * fold, 4);
			% Bootstraps to evaluate the stability of beta.
			bm = [e3 e4];
			d = size(A,1);
			r_by_gene = zeros (size(B,2), cross_repeat_num * fold);
			for p = 1 : cross_repeat_num
				i = datasample(1:d,d,'Replace',false);
				As = A(i,:);
				Bs = B(i,:);
				for q = 0:(fold - 1)
					fprintf('\nRepeat %d - %d', p, q+1);
					idx2 = (mod(1:d,fold) == q);
					Am = As(~idx2,:);
					Bm = Bs(~idx2,:);
					At = As(idx2,:);
					Bt = Bs(idx2,:);
					x = (Am' * Am) \ ( Am' * Bm);

					Bmp = Am * x;
					m_R = mean(arrayfun(@(k) corr(Bmp(k,:)',Bm(k,:)'),1:size(Bm,1),'Uni',1));
					Rm = Bmp - Bm;

					Bp = At * x;
					Rt = Bp - Bt;
					t_R = mean(arrayfun(@(k) corr(Bp(k,:)',Bt(k,:)'),1:size(Bt,1),'Uni',1));
					g_c = arrayfun(@(k) corr(Bp(:,k),Bt(:,k)),1:size(Bt,2),'Uni',1);

					Am = [];
					At = [];

					m_NRSME = sqrt( sum(sum(Rm.^2)) / sum(sum(Bm.^2)) );
					t_NRSME = sqrt( sum(sum(Rt.^2)) / sum(sum(Bt.^2)) );

					Rm = [];
					Rt = [];
					Bm = [];

					idx2 = p * fold + q - fold + 1;
					crv( idx2, 3 ) = m_NRSME;
					crv( idx2, 4 ) = t_NRSME;
					crv( idx2, 1 ) = p;
					crv( idx2, 5 ) = m_R;
					crv( idx2, 6 ) = t_R;
					crv( idx2, 2 ) = q + 1;
					r_by_gene(:,idx2) = g_c';
					
					fprintf('\t%f\t%f\t%f\t%f', m_NRSME, t_NRSME, m_R, t_R);
					
					ep = x(idx);
					bm = [bm ep];
					
				end
			end

			colName ={'Repeat','Iteration','NRMSE_training','NRMSE_test','R_training','R_test'};
			obj.CV_Stat = table(crv(:,1),crv(:,2),crv(:,3),crv(:,4),crv(:,5),crv(:,6),'VariableNames',colName);

			bm = bm(:,3:(cross_repeat_num * fold + 2));
			bm = bm';
			e5 = mean( bm );
			e6 = std(bm);
			e5 = round(e5,5);
			e6 = round(e6,5);
			e5 = e5';
			e6 = e6';
			e7 = round(e5 - e3, 5);
			e8 = round(abs(e7 ./ e3), 4);
			e9 = round(e6 ./ e5, 4);
			colName = {'Gene';'TF';'beta';'beta_pvalue';'CV_beta_mean';'CV_beta_std';'beta_bias';'beta_relative_bias';'relative_std'};
			AllEdges = table(e2,e1,e3,e4,e5,e6,e7,e8,e9,'VariableNames',colName);
			obj.AllEdges = AllEdges;
			obj.Correlation_of_target_gene = r_by_gene;
			fprintf('\n');

		end
	end
end

			
