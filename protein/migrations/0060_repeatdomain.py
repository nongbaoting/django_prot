# Generated by Django 2.2.5 on 2023-12-22 16:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0059_auto_20220814_0039'),
    ]

    operations = [
        migrations.CreateModel(
            name='RepeatDomain',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('dm_names', models.TextField(db_index=True)),
                ('clusters', models.TextField()),
                ('numProtein', models.IntegerField()),
                ('min_len', models.IntegerField()),
                ('max_len', models.IntegerField()),
                ('reprProtein', models.CharField(db_index=True, max_length=255)),
            ],
        ),
    ]
