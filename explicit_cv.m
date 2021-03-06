classdef explicit_cv
	properties
		Training_Sample_Num = [];
		Actual_Target_Exp = [];
		Predicted_Target_Exp = [];
		Test_Sample_Num = [];
		Test_Sample_Stat = [];
		Test_Sample_NRMSE_all = [];
		Test_Sample_Mean_Correlation = [];
	end

	methods (Static)
		function obj=explicit_cv (tf_mtx_training, target_mtx_training, tf_mtx_test, target_mtx_test, test_sample_ids)

			Bm = target_mtx_training;
			Am = [ones(size(tf_mtx_training,1),1) tf_mtx_training];
			obj.Training_Sample_Num = size(Am, 1);

			beta = (Am' * Am) \ (Am' * Bm);
			clear Bm Am tf_mtx_training target_mtx_training; 
			
			B = target_mtx_test;
			A = [ones(size(tf_mtx_test,1),1) tf_mtx_test];
			obj.Test_Sample_Num = size(A, 1);

			id = test_sample_ids;
			id = regexp( test_sample_ids,'([a-zA-Z0-9_\.]+)','once','match');	
			
			Bp = A * beta;
			obj.Predicted_Target_Exp = Bp;
			obj.Actual_Target_Exp = B;
			c = arrayfun(@(k) corr(Bp(k,:)',B(k,:)'),1:size(B,1),'Uni',1);
			obj.Test_Sample_Mean_Correlation = mean(c);

			r = Bp - B;
			rt = r';
			Bt = B';
			fn1 = sqrt(sum(rt.^2));
			fn2 = sqrt(sum(Bt.^2));
			fn3 = fn1 ./ fn2;
			obj.Test_Sample_NRMSE_all = sqrt( (sum(sum(rt.^2))) / (sum(sum(Bt.^2))));

			colName = {'Sample','R_test','NRMSE_test'};
			obj.Test_Sample_Stat = table(id,c',fn3','VariableNames',colName);

		end
	end
end


