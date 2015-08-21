# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0003_technology_motf_file'),
    ]

    operations = [
        migrations.RenameField(
            model_name='technology',
            old_name='motf_file',
            new_name='motif_file',
        ),
    ]
