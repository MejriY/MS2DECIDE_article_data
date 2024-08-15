import json

d_jobs={6:{0:'ID=fb072308be814632bd41594e931aaf7e',
           10:'ID=cc22ff3c58a7407c980d6c87fa6a023e',
           25:'ID=8f5e284c28bc4ca79cdddbfa3b33df68',
           50:'ID=22bef439b05e453a942536e8eae1c34a',
           100:'ID=930b22b654344053ab14436c7d38a273',
           250:'ID=30d7e46783f54ad995edd73ddab36f3a',
           500:'ID=9923c562782140c7bc9abc234d685fda',
           0.02:'ID=ff070334f37044929610491415eca0d3',
          0.1:'ID=ee9aae441e4f4829b06cecdf0f293052'},
        5:{0:'ID=ba85a63134694683967c0a3ffcd378ed',
           10:'ID=d5c3b2ce06a645549fbcc16b1469e577',
           25:'ID=12b769b125f84378b88db33bea4328f3',
           50:'ID=6f9df017c6b84aa2adb29ecf45fa96b4',
           100:'ID=5d71554193ae4f25affbd6bd4af50ec4',
           250:'ID=f6286205af704981b7ef549da2a03d59',
           500:'ID=e1af803246424070a0c0d1b9998680c0',
           0.02:'ID=0d2ef3702412403da022eee9f5ceffe0',
          0.1:'ID=6f2c968c0c20492790eda598bf4821c7'},
        4:{0:'ID=595543cee743484dbc67bb6682a8f73c',
           10:'ID=f26759ac347c46e4b8ce6205da2e53ae',
           25:'ID=3a2b8412d10a4d7abdccfdb6590996c8',
           50:'ID=fb3c140f60f94b50a5db8ead71ec2277',
           100:'ID=d1625ac54fae4f35987e66899b81b8e6',
           250:'ID=5d4fba2cfc7b44938a7a055ae695a944',
           500:'ID=0a7f24fedf604d6199d0495b29ff34da',
           0.02:'ID=1d5bd3667b954f6591748380f4c81d8a',
          0.1:'ID=b2e2569c78b94fb6ad1f47ebe4ce3014'}}
ids_by_first = {k:[v.split('=')[1] for kk,v in d.items()] for k,d in d_jobs.items()}
all_ids = [v for d in ids_by_first.values() for v in d]
print(json.dumps(all_ids))
