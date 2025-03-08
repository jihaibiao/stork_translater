不要强制推送，会把自动发送的action删除  
首先注册github，注册百度翻译高级版，并开通通用文本翻译。注册gmail，开通两步验证，获得应用专用密码，打开iapm。folk我的项目，在setting里面的secrects and variables里面输入各项参数
在本地跑一跑，然后在github的action里面手动跑。

### 期刊数据更新说明
每半年需要手动更新一次 scimagojr.csv 文件，更新步骤：

1. 获取最新数据
   - 访问 https://www.scimagojr.com/journalrank.php
   - 点击 "Download data" 按钮下载最新数据
   - 将下载的文件重命名为 scimagojr.csv

2. 更新仓库
   ```bash
   # 拉取最新代码
   git pull origin main
   
   # 替换文件
   mv 下载目录/scimagojr.csv ./scimagojr.csv
   
   # 提交更新
   git add scimagojr.csv
   git commit -m "更新: scimagojr.csv 数据"
   git push origin main