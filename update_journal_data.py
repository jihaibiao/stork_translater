import pandas as pd
import requests
from pathlib import Path
from datetime import datetime, timedelta
import json

def should_update_data(metadata_file):
    """检查是否需要更新数据"""
    try:
        if not metadata_file.exists():
            return True
            
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
            
        last_update = datetime.fromisoformat(metadata['last_update'])
        # 设置更新间隔为180天（约6个月）
        if datetime.now() - last_update > timedelta(days=180):
            return True
            
        return False
    except Exception:
        return True

def update_scimagojr_data():
    """更新Scimagojr数据"""
    output_file = Path(__file__).parent / 'scimagojr.csv'
    metadata_file = Path(__file__).parent / 'scimagojr_metadata.json'
    
    # 检查是否需要更新
    if not should_update_data(metadata_file) and output_file.exists():
        print("✅ 数据文件仍然有效，无需更新")
        return
    
    try:
        # 检查源文件是否存在
        if not output_file.exists():
            print("❌ 未找到数据源文件 scimagojr.csv")
            return
        
        print("正在处理期刊数据...")
        # 读取CSV数据（注意设置分隔符为分号）
        df = pd.read_csv(output_file, sep=';')
        
        # 只保留需要的列
        df = df[['Title', 'SJR', 'H index', 'SJR Best Quartile']]
        # 重命名列以匹配原代码
        df = df.rename(columns={'SJR Best Quartile': 'Quartile'})
        
        # 保存处理后的数据
        df.to_csv(output_file, index=False)
        
        # 更新元数据
        metadata = {
            'last_update': datetime.now().isoformat(),
            'record_count': len(df)
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f)
        
        print(f"✅ 数据更新成功！保存至: {output_file}")
    except Exception as e:
        print(f"❌ 更新失败: {str(e)}")

if __name__ == "__main__":
    update_scimagojr_data()