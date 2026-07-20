#----------------------------------------------------------------------------------------------------------------
# Data collection
#----------------------------------------------------------------------------------------------------------------

# import libraries
from chembl_webresource_client.settings import Settings
from chembl_webresource_client.new_client import new_client
import pandas as pd
from pathlib import Path

Settings.Instance().TIMEOUT = 5

chembl_id = 'CHEMBL240'

activity = new_client.activity
activities = activity.filter(target_chembl_id=chembl_id).filter(standard_type='IC50')
activities_df = pd.DataFrame(activities)

output_path = Path(__file__).parent.parent.parent / 'data' / 'raw' / 'herg_raw.csv'
activities_df.to_csv(output_path, index=False)