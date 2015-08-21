# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0008_technology_commongraph'),
    ]

    operations = [
        migrations.RenameField(
            model_name='technology',
            old_name='sequence',
            new_name='compared_with',
        ),
    ]
