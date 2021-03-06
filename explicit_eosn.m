classdef explicit_eosn    % eosn: Effect of sample number. To test how the number of training samples affect the predictive power of the predictor model.
	properties
		TestSampleNumber = [];
		stat = [];
	end

	methods (Static)
		function obj=explicit_eosn (x1, x2, TestSampleNum)   % x1: tf gene expression matrix; x2: target gene expression matrix; TestSampleNum: the number of samples held out as test samples.
			
			if nargin < 3
				TestSampleNum = 3000;
			end
			
			obj.TestSampleNumber = TestSampleNum;
				
			B = x2;
			A = [ones(size(x1,1),1) x1];

			i = datasample( 1:size(A,1), TestSampleNum, 'Replace', false);
			idx = ismember( 1:size(A,1), i);

			At = A(idx,:);
			Bt = B(idx,:);

			A = A( ~idx, :);
			B = B( ~idx, :);

			num = floor( size(A, 1)/1000 );
			stat = zeros( 2000, 6);
			SampleNum = 1700;
			i = 0;
			j = 0;

			fprintf('Run\tTrainingSampleNum\tR_training\tNRMSE_training\tR_test  \tNRMSE_test\n');
			while SampleNum <= size(A,1)
				idx = ismember( 1:size(A,1), datasample( 1:size(A,1), SampleNum, 'Replace', false));
				As = A( idx, : );
				Bs = B( idx, : );
				beta = ( As' * As ) \ ( As' * Bs);

				Bp = As * beta;
				R_training = mean ( arrayfun(@(k) corr(Bp(k,:)',Bs(k,:)'),1:size(Bs,1),'Uni',1));
				r = Bp - Bs;
				NRMSE_training = sqrt( (sum(sum(r.^2))) / (sum(sum(Bs.^2))));

				Btp = At * beta;
				R_test = mean( arrayfun(@(k) corr(Btp(k,:)',Bt(k,:)'),1:size(Bt,1),'Uni',1));
				r = Btp - Bt;
				NRMSE_test = sqrt( (sum(sum(r.^2))) / (sum(sum(Bt.^2))));
				
				i = i + 1;
				fprintf(' %d\t            %d\t%f\t%f\t%f\t%f\n',i, SampleNum, R_training, NRMSE_training, R_test, NRMSE_test);
				stat(i,:) = [i SampleNum R_training NRMSE_training R_test NRMSE_test];

				if SampleNum >= 30000
					SampleNum = SampleNum + 5000;
				elseif SampleNum >= 10000
					SampleNum = SampleNum + 2000;
				elseif SampleNum >= 5000
					SampleNum = SampleNum + 1000;
				elseif SampleNum >= 2500
					SampleNum = SampleNum + 500;
				elseif SampleNum >= 2000
					SampleNum = SampleNum + 100;
				elseif SampleNum >= 1800
					SampleNum = SampleNum + 50;
				else
					SampleNum = SampleNum + 25;
				end

				if SampleNum >= size(A,1) && j == 0
					SampleNum = size(A, 1);
					j = 1;
				end

			end

			stat = stat(1:i,:);
			colName = {'Run','TrainingSampleNum','R_training','NRMSE_training','R_test','NRMSE_test'};	
			t = table( stat(:,1), stat(:,2), stat(:,3),stat(:,4),stat(:,5),stat(:,6),'VariableNames',colName);
			obj.stat = t;
		end
	end

end


