# Retro_test

	Python s1_token_process s1_token_process.py /home/kesci/input/competition/train.txt ../data_token train 
	python s1_build_token_vocab.py ../data_token 
	data_gen.sh t2t_token_test1 my_reaction_token ../data_token
	bash data_trainer.sh t2t_token_test1 my_reaction_token
	bash data_decoder.sh t2t_token_test1 my_reaction_token ../data_token test_sources
